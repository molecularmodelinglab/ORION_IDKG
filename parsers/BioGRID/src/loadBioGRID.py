import os
import re
import json
import requests
import pandas as pd
from bs4 import BeautifulSoup
from urllib.parse import urljoin
from zipfile import ZipFile


from Common.utils import GetData, GetDataPullError
from Common.loader_interface import SourceDataLoader
from Common.kgxmodel import kgxnode, kgxedge
from Common.biolink_constants import *


class BioGRIDLoader(SourceDataLoader):
    source_id = 'BIOGRID'
    provenance_id = 'infores:biogrid'
    parsing_version: str = '1.0'

    release_archive_url = "https://downloads.thebiogrid.org/BioGRID/Release-Archive/"
    release_archive_download_url = "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/"
    
    chunk_size = 50000

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode: Run in test mode
        :param source_data_dir: Local directory to store downloaded data
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)

        # ZIP file names that will be downloaded
        self.biogrid_all_zip = 'BIOGRID-ALL-4.4.248.tab3.zip'
        self.biogrid_chemicals_zip = 'BIOGRID-CHEMICALS-4.4.248.chemtab.zip'
        self.zip_files = [self.biogrid_all_zip, self.biogrid_chemicals_zip]
        
        # TXT file names that will be extracted from ZIP files
        self.biogrid_all_file = 'BIOGRID-ALL-4.4.248.tab3.txt'
        self.biogrid_chemicals_file = 'BIOGRID-CHEMICALS-4.4.248.chemtab.txt'
        self.data_files = [self.biogrid_all_file, self.biogrid_chemicals_file]
        self._load_chemical_action_mappings()

    def _load_chemical_action_mappings(self):
        """Load chemical action to predicate mappings from external JSON file"""
        mapping_file = os.path.join(os.path.dirname(__file__), 'chemical_action_mapping.json')
        try:
            with open(mapping_file, 'r') as f:
                self.chemical_action_mapping = json.load(f)
            self.logger.info(f'Loaded {len(self.chemical_action_mapping)} chemical action mappings from {mapping_file}')
        except FileNotFoundError:
            self.logger.error(f'Chemical action mapping file not found: {mapping_file}')
            self.chemical_action_mapping = {}
        except json.JSONDecodeError as e:
            self.logger.error(f'Error parsing chemical action mapping file: {e}')
            self.chemical_action_mapping = {}

    def get_latest_source_version(self) -> str:
        """
        Scrapes the archive directory for the latest release.
        :return: Version string like '4.4.248'
        """
        try:
            response = requests.get(self.release_archive_url)
            soup = BeautifulSoup(response.content, "html.parser")

            version_pattern = re.compile(r'BIOGRID-(\d+\.\d+\.\d+)/')
            versions = []

            for link in soup.find_all("a", href=True):
                href = link["href"]
                match = version_pattern.search(href)
                if match:
                    versions.append(match.group(1))

            if not versions:
                raise Exception("No BIOGRID-x.y.z/ folders found.")

            latest_version = sorted(versions, key=lambda v: tuple(map(int, v.split("."))))[-1]
            return latest_version

        except Exception as e:
            raise GetDataPullError(f"Unable to determine latest version for BioGRID: {e}")

    def get_data(self):
        """
        Downloads and extracts the BIOGRID-ALL and BIOGRID-CHEMICALS ZIP files for the latest version.
        """
        version = self.get_latest_source_version()
        
        # Update file names with the current version
        self.biogrid_all_zip = f'BIOGRID-ALL-{version}.tab3.zip'
        self.biogrid_chemicals_zip = f'BIOGRID-CHEMICALS-{version}.chemtab.zip'
        self.zip_files = [self.biogrid_all_zip, self.biogrid_chemicals_zip]
        
        self.biogrid_all_file = f'BIOGRID-ALL-{version}.tab3.txt'
        self.biogrid_chemicals_file = f'BIOGRID-CHEMICALS-{version}.chemtab.txt'
        
        release_url = urljoin(self.release_archive_download_url, f"BIOGRID-{version}/")
        
        gd: GetData = GetData(self.logger.level)

        for zip_file in self.zip_files:
            full_url = urljoin(release_url, zip_file)
            gd.pull_via_http(full_url, data_dir=self.data_path)
            zip_path = os.path.join(self.data_path, zip_file)
            self._extract_zip_file(zip_path)

        return True
    
    def _extract_zip_file(self, zip_path: str):
        """
        Extracts the contents of a ZIP file to the data directory.
        
        :param zip_path: Path to the ZIP file to extract
        """
        try:
            with ZipFile(zip_path, 'r') as zip_ref:

                zip_ref.extractall(self.data_path)
                self.logger.info(f'Successfully extracted {zip_path}')

                extracted_files = zip_ref.namelist()
                self.logger.info(f'Extracted files: {extracted_files}')
                
        except Exception as e:
            self.logger.error(f'Error extracting ZIP file {zip_path}: {e}')
            raise GetDataPullError(f'Failed to extract ZIP file {zip_path}: {e}')

    def parse_data(self) -> dict:
        """
        Parses the BioGRID data files for graph nodes/edges
        :return: load_metadata
        """
        total_processed = 0
        total_skipped = 0
        
        self.logger.info('Starting BioGRID data parsing...')
        
        # Process protein / gene interactions
        if self.biogrid_all_file:
            protein_file_path = os.path.join(self.data_path, self.biogrid_all_file)
            if os.path.exists(protein_file_path):
                self.logger.info('Processing protein / gene interactions...')
                processed, skipped = self._process_protein_interactions(protein_file_path)
                total_processed += processed
                total_skipped += skipped
                self.logger.info(f'Processed {processed} protein interactions, skipped {skipped}')
        
        # Process chemical interactions
        if self.biogrid_chemicals_file:
            chemical_file_path = os.path.join(self.data_path, self.biogrid_chemicals_file)
            if os.path.exists(chemical_file_path):
                self.logger.info('Processing chemical interactions...')
                processed, skipped = self._process_chemical_interactions(chemical_file_path)
                total_processed += processed
                total_skipped += skipped
                self.logger.info(f'Processed {processed} chemical interactions, skipped {skipped}')
        
        self.logger.info(f'Total processed: {total_processed}, Total skipped: {total_skipped}')
        
        return {
            'num_source_lines': total_processed,
            'unusable_source_lines': total_skipped
        }
    

    def _process_protein_interactions(self, file_path: str) -> tuple:
        """Process protein interactions from BioGRID tab file"""
        processed_count = 0
        skipped_count = 0
        
        try:
            for chunk in pd.read_csv(file_path, sep='\t', chunksize=self.chunk_size, low_memory=False, keep_default_na=False, na_values=['', 'NULL', 'null', 'N/A', 'n/a', 'NA', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A', 'NULL', 'NaN', 'n/a', 'null']):
                for _, row in chunk.iterrows():
                    processed_count += 1
                    
                    try:
                        # Extract data
                        exp_system = row['Experimental System']
                        exp_type = row['Experimental System Type']
                        org_a = row['Organism Name Interactor A']
                        org_b = row['Organism Name Interactor B']
                        
                        # Create CURIEs
                        subject_curie = f"NCBIGene:{row['Entrez Gene Interactor A']}"
                        object_curie = f"NCBIGene:{row['Entrez Gene Interactor B']}"
                        
                        # Determine predicate based on experimental type
                        if exp_type == 'physical':
                            predicate = 'biolink:physically_interacts_with'
                        elif exp_type == 'genetic':
                            predicate = 'biolink:genetically_interacts_with'
                        else:
                            # Log unexpected experimental types for monitoring
                            self.logger.warning(f'Unexpected experimental type: {exp_type} for system: {exp_system}')
                            predicate = 'biolink:interacts_with'
                        
                        # Determine context
                        interaction_context = self._get_interaction_context(org_a, org_b)
                        
                        # Create edge properties
                        edge_props = {
                            'experimental_system': exp_system,
                            'experimental_type': exp_type,
                            'throughput': row['Throughput'],
                            'subject_organism': org_a,
                            'object_organism': org_b,
                            'interaction_context': interaction_context,
                            KNOWLEDGE_LEVEL: KNOWLEDGE_ASSERTION,
                            AGENT_TYPE: MANUAL_AGENT
                        }
                        
                        # Add publication if available
                        if 'PUBMED:' in str(row['Publication Source']):
                            pubmed_id = row['Publication Source'].replace('PUBMED:', '')
                            edge_props[PUBLICATIONS] = [f'PMID:{pubmed_id}']
                        
                        
                        # Create nodes
                        subject_node = kgxnode(
                            identifier=subject_curie,
                            name=row['Official Symbol Interactor A'],
                            categories=[GENE]
                        )
                        self.output_file_writer.write_kgx_node(subject_node)
                        
                        object_node = kgxnode(
                            identifier=object_curie,
                            name=row['Official Symbol Interactor B'],
                            categories=[GENE]
                        )
                        self.output_file_writer.write_kgx_node(object_node)
                        
                        # Create edge
                        edge = kgxedge(
                            subject_id=subject_curie,
                            object_id=object_curie,
                            predicate=predicate,
                            primary_knowledge_source=self.provenance_id,
                            edgeprops=edge_props
                        )
                        self.output_file_writer.write_kgx_edge(edge)
                        
                    except Exception as e:
                        self.logger.warning(f'Error processing protein interaction row {processed_count}: {e}')
                        skipped_count += 1
                        continue
                    
        except Exception as e:
            self.logger.error(f'Error reading protein interactions file: {e}')
            
        return processed_count, skipped_count

    def _generate_chemical_curie(self, row) -> str:
            """Generate chemical CURIE based on source and available identifiers"""
            chemical_source = row.get('Chemical Source')
            chemical_source_id = row.get('Chemical Source ID')
            biogrid_chemical_id = row.get('BioGRID Chemical ID')
            inchi_key = row.get('InChIKey')
            
            # Return BioGRID fallback if no source info available
            if pd.isna(chemical_source) or chemical_source == '-':
                return f"BIOGRID.CHEMICAL:{biogrid_chemical_id}"
            
            # Return BioGRID fallback if no source ID available
            if pd.isna(chemical_source_id) or chemical_source_id == '-':
                return f"BIOGRID.CHEMICAL:{biogrid_chemical_id}"
            
            # Generate CURIE based on chemical source
            source_mapping = {
                'BIOGRID': f"BIOGRID.CHEMICAL:{chemical_source_id}",
                'CHEBI': f"CHEBI:{chemical_source_id}",
                'CHEMBL': f"CHEMBL.COMPOUND:{chemical_source_id}",
                'DRUGBANK': f"DRUGBANK:{chemical_source_id}",
                'PUBCHEM': f"PUBCHEM.COMPOUND:{chemical_source_id}"
            }
            
            if chemical_source in source_mapping:
                return source_mapping[chemical_source]
            
            # Handle CHEMSPIDER and PDB with InChIKey
            if chemical_source in ['CHEMSPIDER', 'PDB']:
                if pd.notna(inchi_key) and inchi_key != '-':
                    return f"INCHIKEY:{inchi_key}"
                else:
                    return f"BIOGRID.CHEMICAL:{biogrid_chemical_id}"
            
            # Unknown source, use BioGRID fallback
            return f"BIOGRID.CHEMICAL:{biogrid_chemical_id}"

    def _process_chemical_interactions(self, file_path: str) -> tuple:
        """Process chemical interactions from BioGRID txt file"""
        processed_count = 0
        skipped_count = 0
        
        try:
            df = pd.read_csv(file_path, sep='\t', low_memory=False, keep_default_na=False, na_values=['', 'NULL', 'null', 'N/A', 'n/a', 'NA', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A', 'NULL', 'NaN', 'n/a', 'null'])
            for _, row in df.iterrows():
                processed_count += 1
                
                try:
                    # Extract basic information
                    organism = row['Organism']
                    action = row['Action']
                    int_type = row['Interaction Type']
                    
                    # Generate CURIEs
                    chemical_curie = self._generate_chemical_curie(row)
                    target_curie = f"NCBIGene:{row['Entrez Gene ID']}"
                    
                    # Map predicate
                    predicate = self.chemical_action_mapping.get(action, 'biolink:chemically_interacts_with')
                    
                    edge_props = {
                        'chemical_action': action,
                        'target_organism': organism,
                        'interaction_type': int_type,
                        'curated_by': row['Curated By'],
                        KNOWLEDGE_LEVEL: KNOWLEDGE_ASSERTION,
                        AGENT_TYPE: MANUAL_AGENT
                    }
    
                    pubmed_id = row.get('Pubmed ID')
                    if pd.notna(pubmed_id) and pubmed_id != '-':
                        edge_props[PUBLICATIONS] = [f"PMID:{pubmed_id}"]
                        
                    atc_codes = row.get('ATC Codes')
                    if pd.notna(atc_codes) and atc_codes != '-':
                        edge_props['atc_codes'] = atc_codes
                    
                    # Create chemical node
                    chemical_node = kgxnode(
                        identifier=chemical_curie,
                        name=row['Chemical Name'],
                        categories=[CHEMICAL_SUBSTANCE]
                    )
                    self.output_file_writer.write_kgx_node(chemical_node)
                    
                    # Create target node
                    target_node = kgxnode(
                        identifier=target_curie,
                        name=row['Official Symbol'],
                        categories=[GENE]
                    )
                    self.output_file_writer.write_kgx_node(target_node)
                    
                    # Create edge
                    edge = kgxedge(
                        subject_id=chemical_curie,
                        object_id=target_curie,
                        predicate=predicate,
                        primary_knowledge_source=self.provenance_id,
                        edgeprops=edge_props
                    )
                    self.output_file_writer.write_kgx_edge(edge)
                        
                except Exception as e:
                    self.logger.warning(f'Error processing chemical interaction row {processed_count}: {e}')
                    skipped_count += 1
                    continue
                    
        except Exception as e:
            self.logger.error(f'Error reading chemical interactions file: {e}')
            raise
        
        return processed_count, skipped_count

    def _get_interaction_context(self, org_a: str, org_b: str) -> str:
        """Determine interaction context"""
        if org_a == org_b:
            return 'within_species'
        else:
            return 'cross_species'
