# Comparative Codon Usage Analysis Pipeline

## Description
This repository contains a Python 3 pipeline to calculate, summarize, and visualize codon usage bias metrics across multiple genomes. 

## Prerequisites
Dependencies can be installed via standard package managers:
`pip install biopython pandas numpy scipy scikit-learn seaborn matplotlib openpyxl`

## Input Data Requirements
1. **File Format:** Sequences must be provided as FASTA files (`.fa`, `.fasta`, or `.fna`).
2. **Directory Structure:** One FASTA file per genome/lineage. The script utilizes the filename as the operational identifier in all output matrices.
3. **Content:** Files must contain Coding DNA Sequences (CDS).
4. **Sequence Quality:** Sequences must be in-frame, beginning with a valid start codon and terminating with a stop codon. The script includes automated QC logging and will flag sequences < 9bp or those lacking valid sense codons.

## Usage
Execute the script via the command line, providing the path to the directory containing your CDS FASTA files:
`python codon_analysis.py /path/to/fasta_folder`

## Outputs
The pipeline processes the input directory and generates the following files:

### 1. Data Matrices & Summaries (.csv, .xlsx)
* **codon_frequencies:** Standard genome-wide codon frequencies.
* **rscu_matrix:** Relative Synonymous Codon Usage (RSCU) values.
* **per_gene_metrics.csv:** Comprehensive master table containing GC content, GC3, ENC, and CAI for every analyzed gene.
* **genome_meta:** Genome-level metadata including total codon counts and mean GC/GC3.
* **codon_summary.csv:** Statistical summaries including variance of mean-RSCU and R-squared values for the regressions.

### 2. Visualizations (.png)
* **codon_freq_heatmap.png & rscu_heatmap.png:** Hierarchically clustered heatmaps (Euclidean distance, average linkage).
* **codon_freq_pca.png & rscu_pca.png:** Principal Component Analysis clustering of genomes based on usage signatures.
* **enc_vs_gc3_plot.png:** Wright's plot visualizing mutational vs. translational selection pressures (colored by genome).
* **gc3_vs_meanrscu_regression.png:** Scatter plots with per-genome linear regression to evaluate the influence of third-position mutational bias.

### 3. Text Reports (.txt)
* **enc_gc3_regression_report.txt:** Exact slope, intercept, and R-squared values for the per-genome ENC vs. GC3 regressions.
* **data_quality_warnings.txt:** QC log identifying any short sequences (<9bp), genes lacking valid codons, or abnormal genome mean RSCUs (only generated if warnings occur).

## Citation References
* **RSCU & CAI:** Sharp, P. M., & Li, W. H. (1987). Nucleic acids research, 15(3), 1281-1295.
* **Methodological note:** CAI is calculated utilizing the complete CDS set of each genome as an internal reference for relative adaptiveness, per Puigbò et al. (2008). Biology Direct, 3, 38.
* **ENC:** Wright, F. (1990). Gene, 87(1), 23-29.