## Functional Enrichment Visualization

### Description
This Python 3 script generates publication-ready bubble plots for ontology and functional enrichment results. It visualizes the significance (log10 P-value) and gene count for enriched terms across multiple genomes simultaneously. The script is pre-configured to output high-resolution figures compliant with standard journal formatting guidelines.

### Prerequisites
`pip install pandas matplotlib Pillow`

### Input Data Requirements
Requires a tab-separated text file (`enrichment_data.txt`) with the following four columns:
1. **Genome:** Population or lineage identifier (must match the color mapping in the script).
2. **Term:** The ontology or functional term description.
3. **Count:** Number of genes associated with the term.
4. **p-value:** The enrichment significance value.

### Usage
Execute the script via the command line:
```bash
python plot_enrichment.py

##Outputs
Generates a multi-lineage bubble plot in three formats:
ontology_vertical_plot.svg
ontology_vertical_plot.pdf
ontology_vertical_plot.tiff