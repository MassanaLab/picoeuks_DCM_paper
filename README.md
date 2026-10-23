# Picoeukaryotic taxa display similar vertical preferences across the Deep Chlorophyll Maximum in contrasting oceanic regions

This repository contains code and data included in:

--------
Béjar-Queijas, S., Obiol, A., Deulofeu-Capo, O., Marín-Vindas, C., Gazulla, C.R., Arrieta, J.M., Ferrera, I., Gasol, J.M., Sebastián, M., and Massana, R. (2026). Picoeukaryotic taxa display similar vertical preferences across the Deep Chlorophyll Maximum in contrasting oceanic regions.

--------

## Scripts to reproduce all figures

1. Download all scripts and data by using:

```
git clone https://github.com/MassanaLab/picoeuks_DCM_paper.git
```

2. Install the following R packages: 

```
speedyseq
ggtree
ape
picante
FSA
phytools
Biostrings
patchwork
tidyverse
broom
e1071
ggh4x
ggnewscale
ggtreeExtra
```

3. Open the [`figures_and_tables.qmd`](scripts/figures_and_tables.qmd) script and run the code to reproduce all figures and tables.

## Data

### Phyloseq objects

Phyloseq objects with taxonomy and metadata can be found in [`data/phyloseq/`](data/phyloseq/).

### Metadata

Other metadata can be found in [`data/metadata/`](data/metadata/)

### Alignments and trees

Fasta files, aligned fasta files and trees can be found in [`data/trees/`](data/trees/).

## Where to find raw data

Raw data can be found here:

|Dataset    |Accession number |
|:----------|:--------------------|
|MIAU       |[PRJEB97768](https://www.ebi.ac.uk/ena/browser/view/PRJEB97768)                  |
|POSEIDON   |[PRJEB97764](https://www.ebi.ac.uk/ena/browser/view/PRJEB97764)                  |
|REMEI      |[PRJEB97766](https://www.ebi.ac.uk/ena/browser/view/PRJEB97766)                  |
