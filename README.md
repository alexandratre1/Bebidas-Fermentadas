# Bebidas-Fermentadas
Python and R workflow used for the metagenomic analysis of fermented beverage samples

# 📖 Overview
This repository contains the bioinformatic workflow used for the metagenomic analysis of Kombucha, Kefir, Pozol and Pulque samples. The project includes:
- Quality control
- Taxonomic classification
- Metagenomic assembly
- Genome binning
- Functional annotation
- Microbial diversity analyses

# 🔍 Repository Structure
```text
fermented-beverages-metagenomics/
├── Source/                      # Software and databases used
├── metadata/                    # Sample metadata
├── environment.yml              # Conda environment file
├── scripts/                     # Python and R scripts for analysis
├── results/                     # Figures and Statistical results
└── README.md                    # Project documentation
```
# 🧬 Metadata
Kombucha
Bio project: PRJNA833075

Kefir
Bio project: PRJNA704713 

Pozol
Bio project: PRJNA648868

Pulque
Bio project: PRJNA603591 

# 🐍 Setup the environment
To recreate the environment with all necessary bioinformatic tools
```bash
# Create the environment
conda env create -f environment.yml

# Activate the environment
conda activate metagenomics
```

# 🧪 Pipeline & Script Descriptions

## 1. Quality Control & Trimming 
- FastQC
- MultiQC
- Trimmomatic
Processes raw Single-End (SE) and Paired-End .fastq.gz reads to remove low-quality bases and adapter sequences using Trimmomatic.
Key parameters: SLIDINGWINDOW:4:20 MINLEN:50

## 2. Taxonomic Classification
Assigns taxonomic reads using Kaiju against the NCBI nr_euk database to identify both bacterial and eukaryotic communities.

## 3. Metagenomic Assembly
Performs de novo assembly on quality-filtered reads using MetaSPAdes to generate contigs for downstream binning.

## 4. Genome Binning
Recovery of MAGs using MaxBin2.

## 4. Downstream Analysis & MAGs
High-quality bins (MAGs) were selected based on 75% completeness

## 5. Taxonomic Classification of MAGs
Taxonomic classification of MAGs using GTDB-Tk.

# 📊 R Analysis and Data Visualization
The R scripts located in scripts_R/ 
Microbial community composition, diversity, and functional profiles of fermented beverage samples were analyzed using R and multiple bioinformatic tools.
Utilize tidyverse and ggplot2 to process taxonomic tables and generate plots.
 If you don't have ggplot2 installed, you can install it in an R session with:
```r
install.packages("ggplot2")
install.packages("tidyverse")
```
The analyses included:
- Taxonomic Analysis
- Alpha Diversity
- Beta Diversity
- Functional Analysis
- Genome Quality Assessment
- Comparative and Visualization Analyses

Example: 

<p align="center">
  <img width="850" alt="Relative Abundance Analysis" src="https://github.com/user-attachments/assets/a91006ae-5685-4abd-8af8-76bac123b1ee" />
</p>
