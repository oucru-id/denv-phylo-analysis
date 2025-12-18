# DENV FHIR Phylogeny Analysis Pipeline

This pipeline processes FHIR bundle JSON files containing pan-serotype Dengue genomics data to generate distance matrices and phylogenetic trees. Full documentation on progress


## Features

- **Direct FHIR Genomics JSON input**
- **Phylogenetic tree (Whole-Genome):** Rectangular, circular, and unrooted. 
- **Transmission network visualization:** Graph and statistical plots (heatmap, violin plot).
- **Clinical metadata integration:** Extracted from Bundle Genomics FHIR.

## Usage

### Requirements

- [Nextflow](https://www.nextflow.io/)
- [mafft](https://github.com/GSLBiotech/mafft)
- iqtree2
- Python 3.8+
- Python packages: `biopython`, `pandas`, `pyvis`, `matplotlib`, `seaborn`

Install Python dependencies:
```bash
pip install biopython pandas networkx pyvis matplotlib seaborn
```

### Run the Pipeline

```bash
nextflow run main.nf
```

### Input

- FHIR Bundle Genomics JSON files in `data/JSON/`
- Reference genome FASTA (DENV1-4)in `data/references/`
