import os
import re
import rdflib
from typing import Tuple, Optional, Dict, List

from Common.utils import GetData, GetDataPullError
from Common.loader_interface import SourceDataLoader
from Common.kgxmodel import kgxnode, kgxedge
from Common.biolink_constants import *


class PathoPhenoDBLoader(SourceDataLoader):
    source_id = 'PATHOPHENODB'
    provenance_id = 'infores:pathophenodb'
    parsing_version: str = '1.0'

    # PathoPhenoDB data URL
    data_url = "http://patho.phenomebrowser.net/media/downloads/patho_pheno_withsymbols.nt"
    
    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode: Run in test mode
        :param source_data_dir: Local directory to store downloaded data
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)
        
        # Data file names
        self.data_file = 'patho_pheno_withsymbols.nt'
        self.cleaned_data_file = 'cleaned_patho_pheno.nt'
        
        # Initialize parser components
        self._init_constants()
        self._init_parser_state()
    
    def _init_constants(self):
        """Initialize RDF and mapping constants"""
        # RDF predicates
        self.HAS_PHENOTYPE = "http://purl.obolibrary.org/obo/RO_0002200"
        self.HAS_PATHOGEN = "http://purl.obolibrary.org/obo/RO_0002556"
        self.HAS_EVIDENCE = "RO_0002558"
        self.TREATED_BY = "http://purl.obolibrary.org/obo/RO_0002302"
        
        # Custom predicates
        self.RESISTANT_DNA = "http://bio2vec.net/RO#resistant_DNAaccession"
        self.RESISTANT_DRUG = "http://bio2vec.net/RO#resistant_to_drug"
        self.ANTIBIOTIC_RESISTANCE = "http://bio2vec.net/RO#antibiotic_resistance"
        self.RESISTANT_PROTEIN = "http://bio2vec.net/RO#resistant_protein"
        
        # Constants
        self.LABEL_SUFFIX = "#label"
        self.ECO_PREFIX = "ECO_"
        self.ASSOCIATION_PREFIX = "bio2vec.net/dis_"
        self.ANNOTATION_RELATION = "SIO_000255"
        
        # Progress reporting
        self.PROGRESS_INTERVAL = 100000
        
        # Predicate mappings to Biolink
        self.predicate_mappings = {
            self.HAS_PHENOTYPE: 'biolink:has_phenotype',
            self.HAS_PATHOGEN: 'biolink:has_infectious_agent',
            self.TREATED_BY: 'biolink:treated_by',
            self.RESISTANT_DRUG: 'biolink:resistant_to',
            self.ANTIBIOTIC_RESISTANCE: 'biolink:associated_with_resistance_to',
        }
        
        # Special predicate mappings with metadata
        self.special_predicate_mappings = {
            self.RESISTANT_DNA: {
                "predicate": 'biolink:resistant_to',
                "extra": {"resistance_target": "DNA"}
            },
            self.RESISTANT_PROTEIN: {
                "predicate": 'biolink:resistant_to', 
                "extra": {"resistance_target": "protein"}
            },
        }
    
    def _init_parser_state(self):
        """Initialize parser state variables"""
        self.assertion_types = {}  # ECO evidence types and their labels
        self.node_labels = {}  # Node URIs and their human-readable labels
        self.node_to_associations = {}  # Maps nodes to their associations
        self.association_to_target_and_type = {}  # Maps associations to targets and relation types
        self.association_to_assertion_type = {}  # Maps associations to evidence types
        self.direct_edges = []  # List of direct edges

    def get_latest_source_version(self) -> str:
        return "Version_3"

    def get_data(self) -> bool:
        """
        Downloads the PathoPhenoDB N-Triples file
        """
        try:
            gd: GetData = GetData(self.logger.level)
            byte_count: int = gd.pull_via_http(self.data_url, self.data_path, saved_file_name=self.data_file)
            
            if not byte_count:
                return False
                
            self.logger.info(f"Downloaded {self.data_file} ({byte_count:,} bytes)")
            return True
            
        except Exception as e:
            raise GetDataPullError(f"Failed to download PathoPhenoDB data: {e}")

    def parse_data(self) -> dict:
        """
        Parses the PathoPhenoDB N-Triples file for graph nodes/edges
        :return: load_metadata
        """
        file_path = os.path.join(self.data_path, self.data_file)
        
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Data file not found: {file_path}")
        
        self.logger.info('Starting PathoPhenoDB data parsing...')
        
        # Parse RDF file (with cleaning if needed)
        rdf_graph = self._parse_rdf_file(file_path)
        if rdf_graph is None:
            raise Exception("Failed to parse RDF file")
        
        self.logger.info(f'Loaded RDF graph with {len(rdf_graph)} triples')
        
        processed_count = 0
        skipped_count = 0
        
        # Process all triples
        for subject, predicate, rdf_object in rdf_graph:
            processed_count += 1
            
            if processed_count % self.PROGRESS_INTERVAL == 0:
                self.logger.info(f"Processed {processed_count:,} triples...")
            
            try:
                subject_str = str(subject)
                predicate_str = str(predicate)
                object_str = str(rdf_object)
                
                if not self._process_rdf_triple(subject_str, predicate_str, object_str):
                    skipped_count += 1
                    
            except Exception as e:
                self.logger.warning(f'Error processing triple {processed_count}: {e}')
                skipped_count += 1
        
        # Build association-based edges
        self._build_association_edges()
        
        # Convert to KGX format and write
        self._write_nodes_and_edges()
        
        self.logger.info(f'Total processed: {processed_count}, Total skipped: {skipped_count}')
        
        return {
            'num_source_lines': processed_count,
            'unusable_source_lines': skipped_count
        }

    def _parse_rdf_file(self, file_path: str) -> Optional[rdflib.Graph]:
        """Parse RDF file and return graph, cleaning if necessary"""
        graph = rdflib.Graph()
        
        try:
            graph.parse(file_path, format="nt")
            return graph
        except Exception as initial_error:
            self.logger.warning(f"Error loading file: {initial_error}")
            self.logger.info("Attempting to clean the file first...")
            
            try:
                cleaned_file_path = self._clean_malformed_ntriples(file_path)
                self.logger.info(f"Loading cleaned file: {cleaned_file_path}")
                graph.parse(cleaned_file_path, format="nt")
                return graph
            except Exception as cleaning_error:
                self.logger.error(f"Failed even after cleaning: {cleaning_error}")
                return None

    def _clean_malformed_ntriples(self, input_file_path: str) -> str:
        """Clean malformed N-Triples and handle multiple IDs per subject/object"""
        output_file_path = os.path.join(self.data_path, self.cleaned_data_file)
        
        self.logger.info(f"Cleaning malformed N-Triples in {input_file_path}")
        
        with open(input_file_path, 'r', encoding='utf-8') as input_file:
            with open(output_file_path, 'w', encoding='utf-8') as output_file:
                
                for line_number, line in enumerate(input_file, 1):
                    
                    if line_number % self.PROGRESS_INTERVAL == 0:
                        self.logger.info(f"  Cleaned {line_number:,} lines...")
                    
                    original_line = line.strip()
                    if not original_line or original_line.startswith('#'):
                        output_file.write(line)
                        continue
                    
                    # Apply common fixes
                    cleaned_line = self._apply_common_fixes(original_line)
                    cleaned_line = self._ensure_proper_ending(cleaned_line)
                    
                    # Try to parse as triple
                    triple_pattern = r'^<([^>]+)>\s+<([^>]+)>\s+<([^>]+)>\s+\.$'
                    triple_match = re.match(triple_pattern, cleaned_line)
                    if not triple_match:
                        output_file.write(cleaned_line + '\n')
                        continue
                    
                    subject, predicate, obj = triple_match.groups()
                    
                    # Handle multiple IDs in subject/object
                    subject_ids = self._split_multiple_ids(subject)
                    object_ids = self._split_multiple_ids(obj)
                    
                    # Write all combinations
                    for subj_id in subject_ids:
                        for obj_id in object_ids:
                            output_file.write(f"<{subj_id}>\t<{predicate}>\t<{obj_id}> .\n")
        
        self.logger.info("Cleaning complete!")
        return output_file_path

    def _apply_common_fixes(self, line: str) -> str:
        """Apply common fixes to malformed lines"""
        fixed_line = line
        fixed_line = re.sub(r'NCBITaxon_\s+(\d+)', r'NCBITaxon_\1', fixed_line)
        fixed_line = re.sub(r'(\w+_)\s+(\d+)', r'\1\2', fixed_line)
        return fixed_line

    def _ensure_proper_ending(self, line: str) -> str:
        """Ensure line has proper N-Triples ending"""
        if not line.endswith(' .'):
            if line.endswith('.'):
                return line[:-1] + ' .'
            else:
                return line + ' .'
        return line

    def _split_multiple_ids(self, uri_string: str) -> List[str]:
        """Split comma-separated URIs and clean them"""
        if ',' in uri_string:
            return [uri.strip() for uri in uri_string.split(',')]
        return [uri_string]

    def _process_rdf_triple(self, subject_uri: str, predicate_uri: str, object_uri: str) -> bool:
        """Process a single RDF triple and route to appropriate handler"""
        
        # Handle label predicates
        if predicate_uri.endswith(self.LABEL_SUFFIX):
            return self._process_label_relation(subject_uri, object_uri)
        
        # Handle annotation relations
        if self.ANNOTATION_RELATION in predicate_uri:
            return self._process_annotation_relation(subject_uri, object_uri)
        
        # Handle association predicates
        if predicate_uri in {self.HAS_PHENOTYPE, self.HAS_PATHOGEN}:
            return self._process_association_relation(subject_uri, predicate_uri, object_uri)
        
        # Handle evidence relations
        if self.HAS_EVIDENCE in predicate_uri:
            return self._process_evidence_relation(subject_uri, object_uri)
        
        # Handle direct edge predicates
        direct_predicates = {
            self.RESISTANT_DNA, self.RESISTANT_DRUG, self.ANTIBIOTIC_RESISTANCE,
            self.RESISTANT_PROTEIN, self.TREATED_BY
        }
        
        if predicate_uri in direct_predicates:
            self._process_direct_relation(subject_uri, predicate_uri, object_uri)
            return True
        
        return True  # Continue processing unknown predicates

    def _process_label_relation(self, subject_uri: str, object_literal: str) -> bool:
        """Process label relations"""
        if self.ECO_PREFIX in subject_uri:
            if subject_uri in self.assertion_types:
                self.logger.warning(f"Duplicate assertion type found: {subject_uri}")
                return False
            self.assertion_types[subject_uri] = object_literal
        else:
            cleaned_label = object_literal.strip('"')
            self.node_labels[subject_uri] = cleaned_label
        return True

    def _process_annotation_relation(self, subject_uri: str, object_uri: str) -> bool:
        """Process annotation relations"""
        if self.ASSOCIATION_PREFIX not in object_uri:
            self.logger.warning(f"Annotation must lead to association: {subject_uri} -> {object_uri}")
            return False
        
        if subject_uri not in self.node_to_associations:
            self.node_to_associations[subject_uri] = []
        
        if object_uri in self.node_to_associations[subject_uri]:
            self.logger.warning(f"Duplicate association: {object_uri} already linked to {subject_uri}")
            return False
        
        self.node_to_associations[subject_uri].append(object_uri)
        return True

    def _process_association_relation(self, subject_uri: str, predicate_uri: str, object_uri: str) -> bool:
        """Process association relations"""
        if self.ASSOCIATION_PREFIX not in subject_uri:
            self.logger.warning(f"Association relation must originate from association: {subject_uri}")
            return False
        
        if subject_uri in self.association_to_target_and_type:
            self.logger.warning(f"Duplicate association entry: {subject_uri}")
            return False
        
        self.association_to_target_and_type[subject_uri] = (object_uri, predicate_uri)
        return True

    def _process_evidence_relation(self, subject_uri: str, object_uri: str) -> bool:
        """Process evidence relations"""
        if self.ASSOCIATION_PREFIX not in subject_uri:
            self.logger.warning(f"Evidence relation must originate from association: {subject_uri}")
            return False
        
        if self.ECO_PREFIX not in object_uri:
            self.logger.warning(f"Evidence relation must target assertion type: {object_uri}")
            return False
        
        self.association_to_assertion_type[subject_uri] = object_uri
        return True

    def _process_direct_relation(self, subject_uri: str, predicate_uri: str, object_uri: str) -> None:
        """Process direct relations"""
        edge_data = {
            'subject': subject_uri,
            'predicate': predicate_uri,
            'object': object_uri,
            'assertion_type': ""
        }
        self.direct_edges.append(edge_data)

    def _build_association_edges(self) -> None:
        """Build edges from processed association data"""
        for node_uri in self.node_to_associations:
            for association_uri in self.node_to_associations[node_uri]:
                if (association_uri in self.association_to_target_and_type and 
                    association_uri in self.association_to_assertion_type):
                    
                    target_uri, relation_type = self.association_to_target_and_type[association_uri]
                    evidence_type = self.association_to_assertion_type[association_uri]
                    
                    edge_data = {
                        'subject': node_uri,
                        'predicate': relation_type,
                        'object': target_uri,
                        'assertion_type': evidence_type
                    }
                    self.direct_edges.append(edge_data)

    def _create_fallback_label(self, uri: str) -> str:
        """Create fallback label from URI when no label is available"""
        return uri.split('/')[-1].replace('_', ' ')

    def _write_nodes_and_edges(self):
        """Convert parsed data to KGX format and write nodes/edges"""
        # Collect all node URIs that appear in edges
        all_node_uris = set()
        for edge_data in self.direct_edges:
            all_node_uris.add(edge_data['subject'])
            all_node_uris.add(edge_data['object'])
        
        # Write all nodes (with labels when available)
        written_nodes = set()
        for node_uri in all_node_uris:
            curie_id, biolink_type = self._convert_uri_to_curie_and_type(node_uri)
            
            if curie_id and biolink_type and curie_id not in written_nodes:
                # Use label if available, otherwise create fallback
                node_label = self.node_labels.get(node_uri, self._create_fallback_label(node_uri))
                
                node = kgxnode(identifier=curie_id, name=node_label, categories=[biolink_type])
                self.output_file_writer.write_kgx_node(node)
                written_nodes.add(curie_id)
        
        # Write edges
        for edge_data in self.direct_edges:
            subject_curie, _ = self._convert_uri_to_curie_and_type(edge_data['subject'])
            object_curie, _ = self._convert_uri_to_curie_and_type(edge_data['object'])
            
            if not subject_curie or not object_curie:
                self.logger.warning(f"Could not convert edge URIs to CURIEs: {edge_data}")
                continue
            
            # Map predicate to Biolink
            biolink_predicate, extra_metadata = self._map_predicate_to_biolink(edge_data['predicate'])
            
            # Create edge properties
            edge_props = {
                KNOWLEDGE_LEVEL: KNOWLEDGE_ASSERTION,
                AGENT_TYPE: MANUAL_AGENT
            }
            
            # Add extra metadata
            if extra_metadata:
                edge_props.update(extra_metadata)
            
            # Add assertion type if available
            if edge_data['assertion_type']:
                assertion_label = self.assertion_types.get(edge_data['assertion_type'], 
                                                         edge_data['assertion_type'])
                edge_props['evidence_type'] = assertion_label
            
            edge = kgxedge(
                subject_id=subject_curie,
                object_id=object_curie,
                predicate=biolink_predicate,
                primary_knowledge_source=self.provenance_id,
                edgeprops=edge_props
            )
            self.output_file_writer.write_kgx_edge(edge)

    def _convert_uri_to_curie_and_type(self, uri: str) -> Tuple[Optional[str], Optional[str]]:
        """Convert URI to CURIE identifier and Biolink type"""
        # Fix common typos
        cleaned_uri = uri.replace("ncbi.nlm.nih.gov/nuccero/", "ncbi.nlm.nih.gov/nuccore/")
        
        # Handle OBO Library URIs
        if "purl.obolibrary.org/obo/" in cleaned_uri:
            suffix = cleaned_uri.split("obo/")[-1]
            curie = suffix.replace("_", ":", 1)
            
            if curie.startswith("DOID:"):
                return curie, DISEASE
            elif curie.startswith("HP:"):
                return curie, PHENOTYPIC_FEATURE
            elif curie.startswith("MP:"):
                return curie, MODEL_ORGANISM_PHENOTYPE
            elif curie.startswith("NCBITaxon:"):
                return curie, ORGANISM_TAXON
            else:
                return curie, ONTOLOGY_CLASS
        
        # Handle NCBI nuccore URIs
        if "ncbi.nlm.nih.gov/nuccore/" in cleaned_uri:
            identifier = cleaned_uri.split("/")[-1]
            return f"GenBank:{identifier}", NUCLEIC_ACID_ENTITY
        
        # Handle NCBI protein URIs
        if "ncbi.nlm.nih.gov/protein/" in cleaned_uri:
            identifier = cleaned_uri.split("/")[-1]
            return f"RefSeq:{identifier}", PROTEIN
        
        # Handle PubChem compound URIs
        if "pubchem.ncbi.nlm.nih.gov/compound/" in cleaned_uri:
            identifier = cleaned_uri.split("/")[-1]
            return f"PUBCHEM.COMPOUND:{identifier}", CHEMICAL_SUBSTANCE
        
        return None, None

    def _map_predicate_to_biolink(self, predicate_uri: str) -> Tuple[str, Dict[str, str]]:
        """Map RDF predicate URI to Biolink predicate and additional metadata"""
        # Check direct mappings
        if predicate_uri in self.predicate_mappings:
            return self.predicate_mappings[predicate_uri], {}
        
        # Check special cases
        if predicate_uri in self.special_predicate_mappings:
            mapping = self.special_predicate_mappings[predicate_uri]
            return mapping["predicate"], mapping["extra"]
        
        # Unknown predicate - log warning and use fallback
        self.logger.warning(f"Could not map predicate to Biolink: {predicate_uri}")
        return RELATED_TO, {}