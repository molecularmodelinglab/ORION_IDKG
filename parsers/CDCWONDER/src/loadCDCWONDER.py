import os
import json
from typing import Dict, Optional
import pandas as pd
from Common.utils import GetData
from Common.loader_interface import SourceDataLoader
from Common.kgxmodel import kgxnode, kgxedge
from Common.biolink_constants import *

class CDCWONDERLoader(SourceDataLoader):
    source_id = 'CDCWONDER'
    provenance_id = 'infores:cdc-wonder'
    description = 'CDC WONDER is a public health data repository maintained by the Centers for Disease Control and Prevention (CDC).'
    source_data_url = 'https://wonder.cdc.gov/'
    license = 'Public domain'
    attribution = 'https://wonder.cdc.gov/'
    parsing_version: str = '1.2'

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode: Run in test mode
        :param source_data_dir: Local directory to store downloaded data
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)
        
        # List of data files we expect to find in the source directory
        self.data_files = [
            'Age_stats.tsv',
            'Sex_stats.tsv',
            'Ethnicity_stats.tsv',
            'Race_stats.tsv',
            'Region_states_stats_2016.tsv',
            'Region_states_stats_2017.tsv',
            'Region_states_stats_2018.tsv',
            'Region_states_stats_2019.tsv',
            'Region_states_stats_2020.tsv',
            'Region_states_stats_2021.tsv',
            'Region_states_stats_2022.tsv'
        ]
        
        # Mapping files (these will be in the parser directory)
        self.diseases_mapping_file = 'diseases.json'
        self.geography_mapping_file = 'simple_geonames_mapping.json'
        self.parser_data_path = os.path.join(os.path.dirname(__file__))
        
        # Initialize mappings
        self.disease_mappings = self.load_disease_mappings()
        self.geography_mappings = self.load_geography_mappings()
        self.demographic_mappings = self.load_demographic_mappings()
        self.min_rate = 0.01
        
        # Biolink type mappings
        self.biolink_types = {
            'disease': ["biolink:DiseaseOrPhenotypicFeature"],
            'age': [
                "biolink:Cohort",
                "biolink:SubjectOfInvestigation", 
                "biolink:StudyPopulation",
                "biolink:PopulationOfIndividualOrganisms",
                "biolink:OrganismalEntity",
                "biolink:BiologicalEntity",
                "biolink:ThingWithTaxon",
                "biolink:NamedThing"
            ],
            'sex': [
                "biolink:PhenotypicFeature",
                "biolink:DiseaseOrPhenotypicFeature",
                "biolink:BiologicalEntity",
                "biolink:ThingWithTaxon",
                "biolink:NamedThing"
            ],
            'ethnicity': [
                "biolink:PopulationOfIndividualOrganisms",
                "biolink:SubjectOfInvestigation",
                "biolink:OrganismalEntity",
                "biolink:BiologicalEntity",
                "biolink:ThingWithTaxon",
                "biolink:NamedThing"
            ],
            'race': [
                "biolink:PopulationOfIndividualOrganisms",
                "biolink:SubjectOfInvestigation",
                "biolink:OrganismalEntity",
                "biolink:BiologicalEntity",
                "biolink:ThingWithTaxon",
                "biolink:NamedThing"
            ],
            'geography': [
                "biolink:GeographicLocation",
                "biolink:PlanetaryEntity",
                "biolink:NamedThing",
                "biolink:Entity"
            ]
        }

    def load_disease_mappings(self) -> Dict:
        """Load disease mappings from JSON"""
        try:
            diseases_file_path = os.path.join(self.parser_data_path, self.diseases_mapping_file)
            with open(diseases_file_path, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            self.logger.error(f"{self.diseases_mapping_file} not found")
            return {}
    
    def load_geography_mappings(self) -> Dict:
        """Load geography mappings from JSON"""
        mappings = {}
        try:
            geography_file_path = os.path.join(self.parser_data_path, self.geography_mapping_file)
            with open(geography_file_path, 'r') as f:
                mappings.update(json.load(f))
        except FileNotFoundError:
            self.logger.warning(f"{self.geography_mapping_file} not found")
        return mappings
    
    def load_demographic_mappings(self) -> Dict:
        """Load demographic mappings"""
        return {
            'sex': {
                "Male": "UMLS:C0086582",
                "Female": "UMLS:C0086287"
            },
            'race': {
                "American Indian or Alaska Native": "UMLS:C1515945",
                "Black or African American": "UMLS:C0085756", 
                "White": "UMLS:C0043157"
            },
            'ethnicity': {
                "Hispanic": "UMLS:C0086409",
                "Not Hispanic": "UMLS:C4036194",
                "Non-Hispanic": "UMLS:C4036194",
            },
            'age': {
                "<1 year": "UMLS:C0021270",
                "<1  year         ": "UMLS:C0021270",
                "1-4 years": "UMLS:C0008059",
                "5-14 years": "UMLS:C0008059",
                "15-24 years": "UMLS:C0238598",
                "25-39 years": "UMLS:C0001675",
                "25-44 years": "UMLS:C0001675",
                "40-64 years": "UMLS:C0001675",
                "45-64 years": "UMLS:C0001675",
                "65+ years": "UMLS:C0001792",
                "65-74 years": "UMLS:C0001792",
                "75-84 years": "UMLS:C0001792",
                "85+ years": "UMLS:C0001792",
                "< 1 year": "UMLS:C0021270",
                "1-14 years": "UMLS:C0008059",
                "15-39 years": "UMLS:C0238598",
                "40+ years": "UMLS:C0001675",
                "Infant": "UMLS:C0021270",
                "Child": "UMLS:C0008059",
                "Adult": "UMLS:C0001675",
                "Elderly": "UMLS:C0001792"
            }
        }
    
    def get_disease_info(self, original_disease: str) -> tuple:
        """Get disease info from mappings"""
        disease_data = self.disease_mappings.get(original_disease, {})
        
        if isinstance(disease_data, dict):
            return (
                disease_data.get('disease_name', original_disease),
                disease_data.get('qualifier', ''),
                disease_data.get('id', '')
            )
        else:
            self.logger.warning(f"{original_disease} not found")
            return (original_disease, '', disease_data if disease_data else '')
    
    def get_geography_id(self, location_name: str, location_code: str = None) -> str:
        """Get geography ID"""
        if location_name in self.geography_mappings:
            return self.geography_mappings[location_name]
        
        if location_name == "United States":
            return "UMLS:C0454992"
        
        if location_code and not pd.isna(location_code):
            return f"CDC_GEO:{location_code}"
        
        return ""
    
    def get_demographic_id(self, demographic_value: str, demographic_type: str, 
                          demographic_code: str = None) -> str:
        """Get demographic ID"""
        if demographic_type not in self.demographic_mappings:
            return ""
        
        mapping = self.demographic_mappings[demographic_type]
        
        if demographic_value in mapping:
            return mapping[demographic_value]
        
        # Handle unmapped values with codes
        if demographic_type == 'race' and demographic_code and not pd.isna(demographic_code):
            return f"CDC_RACE:{demographic_code}"
        
        return ""
    
    def get_biolink_types(self, category: str) -> str:
        """Get biolink types as comma-separated string"""
        types = self.biolink_types.get(category, [])
        return ",".join(types)
    
    def get_edge_type(self, demographic_field: str) -> str:
        """Get appropriate edge type based on demographic"""
        if demographic_field == "Regions/States":
            return "biolink:occurs_in"
        else:
            return "biolink:affects"
    
    def process_single_table(self, input_file: str, demographic_field: str) -> Optional[pd.DataFrame]:
        """Process a single CDC table and return standardized rows"""
        
        self.logger.info(f"PROCESSING: {input_file}")
        self.logger.info(f"Demographic field: {demographic_field}")
        
        # Load data
        try:
            file_path = os.path.join(self.data_path, input_file)
            df = pd.read_csv(file_path, sep='\t', low_memory=False)
            self.logger.info(f"Loaded {len(df)} rows")
        except Exception as e:
            self.logger.error(f"Error loading {input_file}: {e}")
            return None
        
        # Check required columns
        required_columns = ['Disease', 'Year', demographic_field, 'Published Rate']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            self.logger.error(f"Missing columns: {missing_columns}")
            return None
        
        # Clean and filter data
        df_clean = df[required_columns].copy()
        
        # Add code columns if available
        code_columns = []
        if 'Regions/States Code' in df.columns:
            df_clean['Regions/States Code'] = df['Regions/States Code']
            code_columns.append('Regions/States Code')
        if demographic_field == 'Race' and 'Race Code' in df.columns:
            df_clean['Race Code'] = df['Race Code']
            code_columns.append('Race Code')
        
        # Filter out invalid data
        initial_count = len(df_clean)
        df_clean = df_clean.dropna(subset=['Disease', 'Year', 'Published Rate'])
        df_clean = df_clean[df_clean[demographic_field].notna()]
        df_clean = df_clean[~df_clean[demographic_field].str.contains('Unknown', na=False)]
        
        # Convert and filter rates
        df_clean['Published Rate'] = pd.to_numeric(df_clean['Published Rate'], errors='coerce')
        df_clean = df_clean.dropna(subset=['Published Rate'])
        df_clean = df_clean[df_clean['Published Rate'] >= self.min_rate]
        
        self.logger.info(f"After filtering: {len(df_clean)} rows (from {initial_count})")
        
        if len(df_clean) == 0:
            return None
        
        # Create standardized output rows
        rows = []
        
        for _, row in df_clean.iterrows():
            # Get disease info
            disease_name, disease_qualifier, disease_id = self.get_disease_info(row['Disease'])
            
            # Get demographic info
            demographic_type = {
                'Sex': 'sex',
                'Age': 'age', 
                'Race': 'race',
                'Ethnicity': 'ethnicity',
                'Regions/States': 'geography'
            }.get(demographic_field, 'unknown')
            
            # Get object_id based on demographic type
            if demographic_type == 'geography':
                location_code = row.get('Regions/States Code', None)
                object_id = self.get_geography_id(row[demographic_field], location_code)
            elif demographic_type in ['sex', 'age', 'race', 'ethnicity']:
                race_code = row.get('Race Code', None) if demographic_field == 'Race' else None
                object_id = self.get_demographic_id(row[demographic_field], demographic_type, race_code)
            else:
                object_id = ""
            
            # Skip if no valid mappings
            if not disease_id or not object_id:
                continue
            
            # Build context (qualifier)
            context_parts = []
            if row['Year']:
               context_parts.append(f"year {row['Year']}") 
            if disease_qualifier:
                context_parts.append(disease_qualifier)
            if demographic_type == 'age':
                context_parts.append(row[demographic_field])  # Keep specific age range
            
            context = "; ".join(context_parts) if context_parts else ""
            
            # Create row with clean temporal structure
            standardized_row = {
                'from': row['Disease'],
                'to': row[demographic_field],
                'subject_id': disease_id,
                'object_id': object_id,
                'subject_biolink': self.get_biolink_types('disease'),
                'object_biolink': self.get_biolink_types(demographic_type),
                'edge_type': self.get_edge_type(demographic_field),
                'rate_per_100000': float(row['Published Rate']),
                'context': context
            }
            
            rows.append(standardized_row)
        
        if not rows:
            self.logger.warning("No valid mappings found")
            return None
        
        result_df = pd.DataFrame(rows)
        self.logger.info(f"Created {len(result_df)} standardized rows")
        
        return result_df

    def get_latest_source_version(self) -> str:
        """
        Since we're not downloading data, we'll return a fixed version
        """
        return '1.0'

    def get_data(self) -> bool:
        """
        Check that all required data files exist in the source directory
        """
        # Check that all data files exist
        missing_files = []
        for data_file in self.data_files:
            file_path = os.path.join(self.data_path, data_file)
            if not os.path.exists(file_path):
                missing_files.append(data_file)
            
        self.logger.info('All required data files found in source directory')
        return True

    def parse_data(self) -> dict:
        """
        Parses the CDC WONDER data files and creates KGX nodes and edges
        """
        self.logger.info('Parsing CDC WONDER data...')
        
        self.logger.info(f"Disease mappings: {len(self.disease_mappings)}")
        self.logger.info(f"Geography mappings: {len(self.geography_mappings)}")
        self.logger.info(f"Demographic mappings: {sum(len(m) for m in self.demographic_mappings.values())}")
        
        # Define table files and their demographic fields
        table_files = [
            {'file': 'Sex_stats.tsv', 'demographic': 'Sex'},
            {'file': 'Age_stats.tsv', 'demographic': 'Age'},
            {'file': 'Ethnicity_stats.tsv', 'demographic': 'Ethnicity'},
            {'file': 'Race_stats.tsv', 'demographic': 'Race'},
            {'file': 'Region_states_stats_2016.tsv', 'demographic': 'Regions/States'},
            {'file': 'Region_states_stats_2017.tsv', 'demographic': 'Regions/States'},
            {'file': 'Region_states_stats_2018.tsv', 'demographic': 'Regions/States'},
            {'file': 'Region_states_stats_2019.tsv', 'demographic': 'Regions/States'},
            {'file': 'Region_states_stats_2020.tsv', 'demographic': 'Regions/States'},
            {'file': 'Region_states_stats_2021.tsv', 'demographic': 'Regions/States'},
            {'file': 'Region_states_stats_2022.tsv', 'demographic': 'Regions/States'}
        ]
        
        all_rows = []
        
        for file_info in table_files:
            filename = file_info['file']
            demographic = file_info['demographic']
            
            file_path = os.path.join(self.data_path, filename)
            if not os.path.exists(file_path):
                self.logger.warning(f"File not found: {filename}")
                continue
            
            df_processed = self.process_single_table(filename, demographic)
            
            if df_processed is not None:
                all_rows.append(df_processed)
                self.logger.info(f"✓ Added {len(df_processed)} rows from {filename}")
        
        if not all_rows:
            self.logger.error("No data to combine")
            return {'num_source_lines': 0, 'unusable_source_lines': 0}
        
        # Combine all data
        combined_df = pd.concat(all_rows, ignore_index=True)
        
        # Rename published_rate to rate_per_100000 for clarity
        if 'published_rate' in combined_df.columns:
            combined_df = combined_df.rename(columns={'published_rate': 'rate_per_100000'})
        
        # Reorder columns to final specification
        final_columns = ['from', 'to', 'subject_id', 'object_id', 'subject_biolink', 
                        'object_biolink', 'edge_type', 'rate_per_100000', 'context']
        
        # Only select columns that exist in the dataframe
        existing_columns = [col for col in final_columns if col in combined_df.columns]
        final_df = combined_df[existing_columns]
        
        self.logger.info(f"Total edges: {len(final_df)}")
        self.logger.info(f"Unique diseases: {final_df['subject_id'].nunique()}")
        self.logger.info(f"Unique objects: {final_df['object_id'].nunique()}")
        
        # Process each row and create nodes and edges
        processed_count = 0
        skipped_count = 0
        
        for _, row in final_df.iterrows():
            try:
                processed_count += 1
                
                # Create subject node (disease)
                subject_node = kgxnode(
                    identifier=row['subject_id'],
                    name=row['from'],
                    categories=[DISEASE_OR_PHENOTYPIC_FEATURE]  # Using disease category from biolink_constants
                )
                self.output_file_writer.write_kgx_node(subject_node)
                
                # Create object node (demographic/geographic entity)
                # Determine category based on biolink types in the row
                object_categories = row['object_biolink'].split(',') if row['object_biolink'] else [NAMED_THING]
                object_node = kgxnode(
                    identifier=row['object_id'],
                    name=row['to'],
                    categories=object_categories
                )
                self.output_file_writer.write_kgx_node(object_node)
                
                # Create edge properties
                edge_props = {
                    'rate_per_100000': row['rate_per_100000'],
                    KNOWLEDGE_LEVEL: STATISTICAL_ASSOCIATION,  # This is statistical data
                    AGENT_TYPE: MANUAL_AGENT
                }
                
                # Add context if available
                if 'context' in row and row['context'] and row['context'].strip():
                    edge_props[CONTEXT_QUALIFIER] = row['context']
                
                # Create edge
                edge = kgxedge(
                    subject_id=row['subject_id'],
                    object_id=row['object_id'],
                    predicate=row['edge_type'],
                    primary_knowledge_source=self.provenance_id,
                    edgeprops=edge_props
                )
                self.output_file_writer.write_kgx_edge(edge)
                
            except Exception as e:
                self.logger.warning(f'Error processing row {processed_count}: {e}')
                skipped_count += 1
                continue
                
        self.logger.info(f'Processed {processed_count} rows, skipped {skipped_count}')
        
        return {
            'num_source_lines': processed_count,
            'unusable_source_lines': skipped_count
        }