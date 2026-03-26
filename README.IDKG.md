# IDKG — Infectious Disease Knowledge Graph

**Live platform:** https://idkg.mml.unc.edu

## Overview

IDKG is a large-scale, heterogeneous biomedical knowledge graph focused on infectious diseases. It integrates data from 20 curated sources — including SMACC, Heli-SMACC, HVIDB, NPASS, BioGRID, PathoPhenoDB, DrugCentral, ClinicalTrials.gov, and CDC WONDER — into a unified, queryable structure.

The graph is deployed in a **Neo4j** database and follows the **Biolink Model** ontology for node types and edge predicates.

## Statistics

| | |
|---|---|
| Nodes | ~4.6 million |
| Edges | ~26 million |
| Node types | 20 canonical Biolink categories |
| Edge types | 65+ Biolink predicates |

Key node types include Organism Taxon, Small Molecule, Protein, Gene, Disease, Molecular Activity, and Pathway.

## Data Sources

Built using a dedicated fork of the [ORION](https://github.com/RobokopU24/ORION) KG build pipeline:
**https://github.com/molecularmodelinglab/ORION_IDKG**

New parsers were implemented for: SMACC, Heli-SMACC, HVIDB, NPASS, BioGRID, PathoPhenoDB, CDC WONDER, HVPPI.

## Example Query

The following Cypher query explores the molecular context of the SARS-CoV-2 replicase polyprotein R1AB — finding molecular activities it participates in, its connection to Remdesivir, and any COVID-19 associations:

```cypher
MATCH (n0_0:`biolink:MolecularActivity`)-[r0_0:`biolink:has_input`|`biolink:has_output`]-(n1_0:`biolink:GeneOrGeneProduct`)
WHERE n1_0.name_lower IN ['r1ab_sars2 replicase polyprotein 1ab (sprot)']

OPTIONAL MATCH (rem:`biolink:ChemicalOrDrugOrTreatment`)-[r_affects:`biolink:affects`]-(n1_0)
WHERE rem.name_lower IN ['remdesivir']

OPTIONAL MATCH (covid:`biolink:DiseaseOrPhenotypicFeature`)-[r_part:`biolink:has_part`]-(n1_0)
WHERE covid.name_lower IN ['covid-19']

RETURN n0_0, r0_0, n1_0, rem, r_affects, covid, r_part
```

This query returns the subgraph connecting:
- **R1AB** (SARS-CoV-2 replicase polyprotein) — a key viral protein
- **Molecular activities** in which R1AB serves as input or output
- **Remdesivir** — an antiviral drug affecting R1AB
- **COVID-19** — disease associated with R1AB via `has_part`


## Contact

Molecular Modeling Lab, UNC Chapel Hill — https://molecularmodelinglab.github.io
