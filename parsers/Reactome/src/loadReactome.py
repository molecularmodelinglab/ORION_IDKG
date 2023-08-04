import os
import argparse
import requests
import re
import json
import neo4j

from bs4 import BeautifulSoup
from Common.loader_interface import SourceDataLoader
from Common.kgxmodel import kgxnode, kgxedge
from Common.neo4j_tools import Neo4jTools
from Common.prefixes import REACTOME, NCBITAXON, GTOPDB, UNIPROTKB, CHEBI, KEGG_COMPOUND, KEGG_GLYCAN, PUBCHEM_COMPOUND, NCBIGENE, CLINVAR
from Common.utils import GetData


SUBJECT_COLUMN = 0
PREDICATE_COLUMN = 1
OBJECT_COLUMN = 2
INCLUDE_COLUMN = 3


PREDICATE_PROP_COLUMN = -1

PREDICATE_MAPPING = {"compartment": "biolink:occurs_in",
                     "output": "biolink:has_output",
                     "input": "biolink:has_input",
                     "hasEvent": "biolink:contains_process",
                     "precedingEvent": "biolink:precedes",
                     "activeUnit": "biolink:actively_involves",
                     "hasComponent": "biolink:has_part",
                     "catalystActivity": "biolink:actively_involves",
                     "cellType": "biolink:located_in",
                     "goBiologicalProcess": "biolink:subclass_of",
                     "disease": "biolink:disease_has_basis_in"}


# TODO - use something like this instead of manipulating strings for individual cases
#  reactome databaseName: normalizer preferred curie prefix
CURIE_PREFIX_MAPPING = {
    'UniProt': UNIPROTKB,
    'Guide to Pharmacology': GTOPDB,
    'ChEBI': CHEBI,
    'REACT': REACTOME,
    'COMPOUND': KEGG_COMPOUND,
    'PubChem Compound': PUBCHEM_COMPOUND,
    'PubChem Substance': PUBCHEM_COMPOUND,
    'KEGG Glycan': KEGG_GLYCAN,
    'NCBI Entrez Gene': NCBIGENE,
    'ClinVar': CLINVAR

}

#TODO Make SUPPOSE_PROPERTIES into a place to write nodes that should recieve extra properties.
#Then create different ones for different types, like Summation, CatalystActivity, etc.
SUPPOSE_PROPERTIES = ("Summation",)

# Normalized node are normalized once formatted as curie
# ReactionLikeEvent: ("Pathway", "Event", "BlackboxEvent", "FailedReaction", "Depolymerisation", "Polymerisation")
NORMALIZED_NODES = ("ReactionLikeEvent","Pathway", "Event", "BlackboxEvent", "FailedReaction", "Depolymerisation", "Polymerisation")

# These nodes have no reactome stID, however their ID is the combination of two attributes on them
# GO_Term: "Compartment", "GO_BiologicalProcess", "GO_MolecularFunction", "GO_CellularComponent",
# ExternalOntology: "Disease", "CellType"
ON_NODE_MAPPING = ("GO_Term", "Species", "ExternalOntology","ReferenceTherapeutic","ReferenceMolecule","ReferenceSequence","Disease")#databaseName+:+identifier

# Cross mapping nodes require looking for a second node with a normalizable database identifier.
# ReferenceIsoform is the reference ID for EntityWITHACCESSIONEDSEQUENCE(PROTEIN)
# GenomeEncodedEntity: EntityWITHACCESSIONEDSEQUENCE(Protein), EntityWITHACCESSIONEDSEQUENCE(Gene and transcript), EntityWITHACCESSIONEDSEQUENCE(DNA), EntityWITHACCESSIONEDSEQUENCE(RNA)
CROSS_MAPPING = ('EntityWithAccessionedSequence','GenomeEncodedEntity', 'SimpleEntity', 'Drug', 'Complex', 'Polymer')

TO_WRITE = ('Provenance/Include', 'Attribute/Include') #Descriptive features of other existing nodes eg Summation
TO_MAP = ('IDMapping/Include', ) # Maps the external identifier of a node from another node
TO_INCLUDE = ('Include',)
RDF_EDGES_TO_INCLUDE = ('RDF_edges/Include',)
TO_SWITCH = ('Include/SwitchSO', )

##############
# Class: Reactome loader
#
# By: Ola Olasunkanmi, Jon-Michael Beasley, and Evan Morris
# Date: 7/20/2023
# Desc: Class that loads/parses the reactome data.
##############
class ReactomeLoader(SourceDataLoader):

    # Setting the class level variables for the source ID and provenance
    source_id: str = 'Reactome'
    provenance_id: str = 'infores:reactome'
    description = "Reactome is a free, open-source, curated and peer-reviewed pathway database"
    source_data_url = "https://reactome.org/"
    license = "https://reactome.org/license"
    attribution = "https://academic.oup.com/nar/article/50/D1/D687/6426058?login=false"
    parsing_version = '1.1'

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode - sets the run into test mode
        :param source_data_dir - the specific storage directory to save files in
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)
        self.version_url: str = 'https://reactome.org/about/news'

        self.neo4j_dump_file = 'reactome.graphdb.dump'
        self.data_url = 'https://reactome.org/download/current/'
        self.data_files = [self.neo4j_dump_file]

        self.triple_file: str = 'reactomeContents_CriticalTriples.csv'
        self.triple_path = os.path.dirname(os.path.abspath(__file__))

        self.dbid_to_node_id_lookup = {}

    def get_latest_source_version(self) -> str:
        """
        gets the latest available version of the data
        :return:
        """
        # load the web page for CTD
        html_page: requests.Response = requests.get(self.version_url)
        if html_page.status_code == 200:
            # get the html into a parsable object
            resp: BeautifulSoup = BeautifulSoup(html_page.content, 'html.parser')

            # find the version tag
            a_tag: BeautifulSoup.Tag = resp.find_all()
            for tag in a_tag:
                match = re.search(r'V\d+', tag.get_text())
                if match:
                    version = f'{match.group(0)!r}'
                    return version.strip("'")
            return 'version_broken'
        else:
            html_page.raise_for_status()

    def get_data(self) -> bool:
        gd: GetData = GetData(self.logger.level)
        for dt_file in self.data_files:
            gd.pull_via_http(f'{self.data_url}{dt_file}',
                             self.data_path)
        return True
 
    def parse_data(self):
        neo4j_tools = Neo4jTools()
        neo4j_tools.load_backup_dump(f'{self.data_path}/{self.neo4j_dump_file}')
        neo4j_tools.start_neo4j()
        neo4j_tools.wait_for_neo4j_initialization()
        neo4j_driver = neo4j_tools.neo4j_driver

        parse_metadata = self.extract_data(neo4j_driver=neo4j_driver)

        neo4j_driver.close()
        neo4j_tools.stop_neo4j()
        return parse_metadata

    def extract_data(self, neo4j_driver) -> dict:
        #These Mapping files contains the Ids that normalizes, so instead of the triples alone, we include the mapping too
        queries_to_include = []
        queries_to_map = []
        queries_to_write = []

        with open(os.path.join(self.triple_path, self.triple_file), 'r') as inf:
            lines = inf.readlines()
            lines_to_process = lines[1:]
        
        for line in lines_to_process:
            line = line.strip().split(',')
            if line[INCLUDE_COLUMN] in RDF_EDGES_TO_INCLUDE:
                queries_to_include.append(
                    self.rdf_edge_mapping(line[SUBJECT_COLUMN], line[PREDICATE_COLUMN], line[OBJECT_COLUMN]))
            elif line[INCLUDE_COLUMN] in TO_SWITCH:
                cypher_query = self.map_ids(line[SUBJECT_COLUMN], line[PREDICATE_COLUMN], line[OBJECT_COLUMN], switch=True)
                queries_to_include.append(cypher_query)
            elif line[INCLUDE_COLUMN] in TO_INCLUDE:
                cypher_query = self.map_ids(line[SUBJECT_COLUMN], line[PREDICATE_COLUMN], line[OBJECT_COLUMN], switch=False)
                queries_to_include.append(cypher_query)
            elif line[INCLUDE_COLUMN] in TO_MAP:
                queries_to_map.append(self.map_ids(line[SUBJECT_COLUMN], line[PREDICATE_COLUMN], line[OBJECT_COLUMN]))
            else:
                if line[INCLUDE_COLUMN] in TO_WRITE:
                    pass
                    #TODO Fix the write properties function used in the line below to write properties like support publications and summations.
                    #queries_to_write.append(self.write_properties(line[SUBJECT_COLUMN], line[PREDICATE_COLUMN], line[OBJECT_COLUMN]))

        self.logger.info(f'mapped queries: {json.dumps(queries_to_map, indent=4)}')
        self.logger.info(f'include queries: {json.dumps(queries_to_include, indent=4)}')
        self.logger.info(f'write queries: {json.dumps(queries_to_write, indent=4)}')

        record_counter: int = 0
        skipped_record_counter: int = 0
        with neo4j_driver.session() as session:

            # TODO Suggested implementation
            # reference_entity_lookup = self.get_reference_entity_mapping(neo4j_session=session)
            # self.dbid_to_node_id_lookup.update(reference_entity_lookup)

            #Map the id from the nodes databaseName+:+identifier
            # for node in ALL_ON_NODE_MAPPING:
            #     node_mapping_cypher = self.on_node_mapping(node)
            #     self.logger.info(f'Running query {node_mapping_cypher}')
            #     session.run(node_mapping_cypher)

            # for node in CROSS_MAPPING:
            #     cross_mapping_cypher = self.cross_map_ids(node)
            #     self.logger.info(f'Running query {cross_mapping_cypher}')
            #     session.run(cross_mapping_cypher)

            #Map other properties from the nodes summation catalystactivity ...
            # PS: Not currently in use since the content file contains only 'include'
            for cypher_query in queries_to_write:
                if cypher_query:
                    session.run(cypher_query)

            #Map the id from the other nodes respectively
            # PS: Not currently in use since the content file contains only 'include'
            for cypher_query in queries_to_map:
                if cypher_query:
                    session.run(cypher_query)

            #Finally, run the cypher to get results - nodes and edges
            for cypher_query in queries_to_include:
                self.logger.info(f'Running query ({cypher_query})...')
                result: neo4j.Result = session.run(cypher_query)
                self.logger.info(f'Cypher query ({cypher_query}) complete.')
                # self.logger.info(f'Sample results: {result[:5]}')
                record_count, skipped_record_count = self.write_neo4j_result_to_file(result)
                record_counter += record_count
                skipped_record_counter += skipped_record_count

        parse_metadata = {
            'num_source_lines': record_counter,
            'unusable_source_lines': skipped_record_counter
        }
        return parse_metadata

    # TODO Suggested implementation
    #  this has not been tested but I would make a cypher call like this to find all of the reference entity mappings
    """
    def get_reference_entity_mapping(self, neo4j_session):
        reference_entity_mapping = {}
        reference_entity_query = 'MATCH (a)-[:referenceEntity|crossReference]->(b) return a,b'
        reference_entity_result = neo4j_session.run(reference_entity_query)
        for record in reference_entity_result:
            # do all reference entities have databaseName and identifier? if so it's easy -
            record_data = record.data()
            curie_prefix = CURIE_PREFIX_MAPPING[record_data['b']['databaseName']]
            reference_entity_mapping[record_data['a']['identity']] = f"{curie_prefix}:{record_data['b']['identifier']}"
        return reference_entity_mapping
    """
    def write_neo4j_result_to_file(self, result: neo4j.Result):
        record_count = 0
        skipped_record_count = 0
        for record in result:
            record_data = record.data()

            node_a_id = self.process_node_from_neo4j(record_data['a_id'], record_data['a'], record_data['a_labels'])
            node_b_id = self.process_node_from_neo4j(record_data['b_id'], record_data['b'], record_data['b_labels'])
            if node_a_id and node_b_id:
                self.process_edge_from_neo4j(node_a_id, record_data['r_type'], node_b_id)
                record_count += 1
            else:
                skipped_record_count += 1
        return record_count, skipped_record_count
    
    def process_node_from_neo4j(self, node_identity, node: dict, node_labels: list = None):
        #self.logger.info(f'processing node: {node_identity}')
        node_id = None
        if any(x in ON_NODE_MAPPING for x in node_labels) or any(x in CROSS_MAPPING for x in node_labels):
            #On-node/Same node ID Mapping eg GO, Disease ...
            if "Species" in node_labels:
                node_id = f"{NCBITAXON}:{node['taxId']}"
            elif "GO_Term" in node_labels:
                try:
                    node_id = f"{CURIE_PREFIX_MAPPING[node['databaseName']]}:{node['accession']}"
                except: #Allow the mapping to fail so we can see what isn't normalizing in failure logs.
                    node_id = f"{node['databaseName']}:{node['accession']}"
            else:
                try:
                    node_id = f"{CURIE_PREFIX_MAPPING[node['databaseName']]}:{node['identifier']}"
                except: #Allow the mapping to fail so we can see what isn't normalizing in failure logs.
                    node_id = f"{node['databaseName']}:{node['identifier']}"
        elif any(x in NORMALIZED_NODES for x in node_labels): 
            node_id = f"{REACTOME}:{node['stId']}"
        if not node_id:
            self.logger.warning(f'A node ID could not be mapped for: {node_identity} (ref_map_labels: {node_labels})')
            return None
        
        node_properties = {}
        if 'definition' in node:
            node_properties['definition'] = node['definition']
        if 'url' in node:
            node_properties['url'] = node['url']
        node_name = node['displayName'] if 'displayName' in node else ''
        node_to_write = kgxnode(node_id, name=node_name, nodeprops=node_properties)
        self.output_file_writer.write_kgx_node(node_to_write)
        return node_id

        # TODO we should replace the previous line with a consolidated node identifier mapping section, and remove
        #  the in-neo4j mapping cypher calls. This should follow a hierarchy of preferred identifier mappings,
        #  based on which will normalize the best etc, something like this:
        """
        # did we map this node already?
        if node['dbId'] in self.dbid_to_node_id_lookup:
            node_id = self.dbid_to_node_id_lookup[node['dbId']]
        else:
            # if the node has a databaseName and identifier can we map it easily? - which cases does this not work for?
            if 'databaseName' in node:
                curie_prefix = CURIE_PREFIX_MAPPING.get(node['databaseName'], None)
                if curie_prefix is None:
                    self.logger.warning(f'Could not find a curie prefix mapping for databaseName {node["databaseName"]}')
                else:
                    node_id = f'{curie_prefix}{node["identifier"]}'
            # if no databaseName mapping see if we can use a different way?
            if node_id is None:
                if node['dbId'] in self.reference_entity_lookup:   # from the suggested implementation above
                    node_id = self.reference_entity_lookup[node['dbId']]    
                elif some other identifier or condition:
                    # for example do we need to use stdId?
                    node_id = xxxxxxxxx
                elif GO Term mapping:
                    node_id = xxxxxxxxx
                else if node_labels:
                    if 'Species' in node_labels:
                        node_id = f'{NCBITAXON}:{node["taxId"]}':
                    elif other node type based mappings?
                else:
                    self.logger.warning('uh oh - couldn't find an identifier for node {node}')
                    return None
        # add whatever we found to the lookup map
        self.dbid_to_node_id_lookup[node['dbId']] = node_id
        """

    def process_edge_from_neo4j(self, subject_id: str, relationship_type: str, object_id: str):
        predicate = PREDICATE_MAPPING.get(relationship_type, None)
        if predicate:
            output_edge = kgxedge(
                subject_id=subject_id,
                object_id=object_id,
                predicate=predicate,
                primary_knowledge_source=self.provenance_id
            )
            self.output_file_writer.write_kgx_edge(output_edge)
        else:
            self.logger.warning(f'A predicate could not be mapped for relationship type {relationship_type}')

    def write_properties(self, s, p, o):
        # TwoWords -> twoWords
        def repl_func(match):
            return match.group(1).lower()
        cypher = None
        #Cross Property Mapping eg Summation
        #TODO Need to work on this function to correctly get properties for certain nodes, like Summations or reference publications.
        if (s in SUPPOSE_PROPERTIES) & (o in SUPPOSE_PROPERTIES):
            if (o == "Summation") & (s == "Summation"):
                cypher = f"MATCH (props)<-[r0:summation]-(a:{s})-[r1:{p}]->(b:{o})-[r2:summation]-" \
                        f"RETURN ref as a, labels(a) as a_labels, id(a) as a_id, type(r1) as r_type, b, labels(b) as b_labels, id(b) as b_id"

        elif o in SUPPOSE_PROPERTIES:
            if o == "Summation":
                cypher = f"MATCH (a:{s})-[r1:{p}]->(b:{o})-[r0:summation]->(props) " \
                        f"RETURN ref as a, labels(a) as a_labels, id(a) as a_id, type(r1) as r_type, b, labels(b) as b_labels, id(b) as b_id"

        elif s in SUPPOSE_PROPERTIES:
            if s == "Summation":
                cypher = f"MATCH (props)<-[r0:summation]-(a:{s})-[r1:{p}]->(b:{o}) " \
                        f"RETURN ref as a, labels(a) as a_labels, id(a) as a_id, type(r1) as r_type, b, labels(b) as b_labels, id(b) as b_id"
        
        return cypher
    
    #TODO Either delete the following function or use it in process_nodes_from_neo4j to do on-node mapping.
    """
    def on_node_mapping(self, node_identity, node, node_labels):
        #On-node/Same node ID Mapping eg GO, Disease ...
        #self.logger.info(f'processing node: {node_identity}')
        node_id = None
        if any(x in ON_NODE_MAPPING for x in node_labels) or any(x in CROSS_MAPPING for x in node_labels):
            #On-node/Same node ID Mapping eg GO, Disease ...
            if "Species" in node_labels:
                node_id = f"{NCBITAXON}:{node['taxId']}"
            elif "GO_Term" in node_labels:
                try:
                    node_id = f"{CURIE_PREFIX_MAPPING[node['databaseName']]}:{node['accession']}"
                except: #Allow the mapping to fail so we can see what isn't normalizing in failure logs.
                    node_id = f"{node['databaseName']}:{node['accession']}"
            else:
                try:
                    node_id = f"{CURIE_PREFIX_MAPPING[node['databaseName']]}:{node['identifier']}"
                except: #Allow the mapping to fail so we can see what isn't normalizing in failure logs.
                    node_id = f"{node['databaseName']}:{node['identifier']}"
        elif any(x in NORMALIZED_NODES for x in node_labels): 
            node_id = f"{REACTOME}:{node['stId']}"
        if not node_id:
            self.logger.warning(f'A node ID could not be mapped for: {node_identity} (ref_map_labels: {node_labels})')
            return None
    """
    def map_ids(self, s, p, o, switch=False):
        # ReferenceID Mappings eg EntityWithAccessionedSequence ->ReferenceSequence ......
        cypher = None
        # OBJECT
        if (s in CROSS_MAPPING) & (o in CROSS_MAPPING):
            if switch == True:
                cypher = f"MATCH (refa)-[r0:referenceEntity|crossReference]-(b:{s})-[r1:{p}]->(a:{o})-[r0:referenceEntity|crossReference]-(refb) " \
                        f"RETURN refa as a, labels(a) as a_labels, id(a) as a_id, type(r1) as r_type, refb as b, labels(b) as b_labels, id(b) as b_id"
            else:
                cypher = f"MATCH (refa)-[r0:referenceEntity|crossReference]-(a:{s})-[r1:{p}]->(b:{o})-[r0:referenceEntity|crossReference]-(refb) " \
                        f"RETURN refa as a, labels(a) as a_labels, id(a) as a_id, type(r1) as r_type, refb as b, labels(b) as b_labels, id(b) as b_id"
        elif s in CROSS_MAPPING:
            if switch == True:
                cypher = f"MATCH (ref)-[r0:referenceEntity|crossReference]-(b:{s})-[r1:{p}]->(a:{o}) " \
                        f"RETURN ref as a, labels(a) as a_labels, id(a) as a_id, type(r1) as r_type, b, labels(b) as b_labels, id(b) as b_id"
            else:
                cypher = f"MATCH (ref)-[r0:referenceEntity|crossReference]-(a:{s})-[r1:{p}]->(b:{o}) " \
                        f"RETURN ref as a, labels(a) as a_labels, id(a) as a_id, type(r1) as r_type, b, labels(b) as b_labels, id(b) as b_id"
        elif o in CROSS_MAPPING:
            if switch == True:
                cypher = f"MATCH (b:{s})-[r1:{p}]->(a:{o})-[r0:referenceEntity|crossReference]-(ref) " \
                        f"RETURN a, labels(a) as a_labels, id(a) as a_id, type(r1) as r_type, ref as b, labels(b) as b_labels, id(b) as b_id"
            else:
                cypher = f"MATCH (a:{s})-[r1:{p}]->(b:{o})-[r0:referenceEntity|crossReference]-(ref) " \
                        f"RETURN a, labels(a) as a_labels, id(a) as a_id, type(r1) as r_type, ref as b, labels(b) as b_labels, id(b) as b_id"
        if not cypher:
            # None of the s or object has external mapping file
            if switch == True:
                cypher = f"MATCH (b:{s})-[r:{p}]->(a:{o}) " \
                        f"RETURN a, labels(a) as a_labels, id(a) as a_id, type(r) as r_type, b, labels(b) as b_labels, id(b) as b_id"
            else:
                cypher = f"MATCH (a:{s})-[r:{p}]->(b:{o}) " \
                        f"RETURN a, labels(a) as a_labels, id(a) as a_id, type(r) as r_type, b, labels(b) as b_labels, id(b) as b_id"
        return cypher

    def rdf_edge_mapping(self, s, p, o):
        cypher = None
        # For any triple like
            # CatalystActivity-[activity]->Go_MolecularFunction and,
            # Reaction-[catalystActivity]->CatalystActivity
                # Collapse the 2 hops and get: Reaction-[activity]-Go_MolecularFunction
        # CatalystActivity to catalystActivity
        if s =='CatalystActivity':
            cypher= f"MATCH (a)-[:catalystActivity]->(x:{s})-[r:activity]->(b:{o}) " \
                    f"RETURN a, labels(a) as a_labels, id(a) as a_id, type(r) as r_type, b, labels(b) as b_labels, id(b) as b_id"
        elif o =='CatalystActivity':
            cypher = f"MATCH (a:{s})-[r:{p}]-(x:{o})-[rx:catalystActivity]-(b) return a, r, b" \
                     f"RETURN a, labels(a) as a_labels, id(a) as a_id, type(r) as r_type, b, labels(b) as b_labels, id(b) as b_id"

        return cypher

if __name__ == '__main__':
        # create a command line parser
    ap = argparse.ArgumentParser(description='Load Reactome data files and create KGX import files.')

    ap.add_argument('-r', '--data_dir', required=True, help='The location of the KGs data file')

    # # parse the arguments
    args = vars(ap.parse_args())

    # # this is the base directory for data files and the resultant KGX files.
    rct_data_dir: str = args['data_dir']

    # get a reference to the processor
    ch_ldr = ReactomeLoader()

    # load the data files and create KGX output
    ch_ldr.load(rct_data_dir + '/nodes.jsonl', rct_data_dir + '/edges.jsonl')
