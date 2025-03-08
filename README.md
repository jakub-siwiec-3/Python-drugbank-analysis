# DrugBank Data Analysis

## Description
This project analyzes data from an XML file provided by DrugBank. It processes information about drugs, synonyms, products, pathways, targets, and drug interactions. The results are presented as data tables (using pandas DataFrames) and visualizations (graphs and charts).

## Requirements
- Python 3.7+
- Libraries:
  - `pandas`
  - `networkx`
  - `matplotlib`
  - `xml.etree.ElementTree` (standard library)
  - `random` (standard library)

## Function Descriptions:
- extract_drug_info – Parses basic drug information (ID, name, description, state, food interactions).
- extract_synonyms – Retrieves synonyms for each drug.
- draw_synonym_graph – Draws a complete graph connecting a drug and its synonyms.
- extract_products – Extracts product information associated with each drug.
- extract_pathways – Extracts pathway information associated with the drugs.
- draw_pathway_bipartite_graph – Draws a bipartite graph showing associations between pathways and drugs.
- plot_drug_histogram – Generates a histogram showing the frequency of drug occurrences in different pathways.
- extract_targets – Extracts information about the drug targets.
- plot_targets_distribution – Creates a pie chart of target distribution by cellular location.
- extract_drug_groups and plot_drug_groups – Extracts drug groups (approved, experimental, etc.) and presents the results using pie charts.
- extract_targets_details and extract_drug_interactions – Retrieves detailed information about targets and drug interactions.
- build_gene_graph, extract_gene_info, and plot_gene_nucleotide_structure – For a specified gene (e.g., IFNAR2), builds a graph linking drugs and products, collects basic gene information, and visualizes the nucleotide sequence.
- extract_classification and plot_classification – Parses drug classification data and creates bar charts showing the number of drugs in different categories.
