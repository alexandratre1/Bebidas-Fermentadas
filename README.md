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

 ```bash
fastqc *_trimmed.fastq.

%%%	Unzip documents

for filename in *.zip
do
unzip $filename
done
```

- MultiQC

 ```bash
multiqc
```

- Trimmomatic
```bash
scripts/trimming_SE
scripts/trimming_PE
```
Processes raw Single-End (SE) and Paired-End .fastq.gz reads to remove low-quality bases and adapter sequences using Trimmomatic.
Key parameters: SLIDINGWINDOW:4:20 MINLEN:50

## 2. Taxonomic Classification
```bash
scripts/taxonomic_analysis
```
Assigns taxonomic reads using Kaiju against the NCBI nr_euk database to identify both bacterial and eukaryotic communities.

## 3. Metagenomic Assembly
```bash
scripts/mags_assembly
```
Performs de novo assembly on quality-filtered reads using MetaSPAdes to generate contigs for downstream binning.

## 4. Genome Binning
Recovery of MAGs using MaxBin2.

```bash
run_MaxBin.pl \
-contig assembly/final.contigs.fa \
-abund abundance.txt \
-out maxbin2/bin \
-thread 16
```

## 4. Downstream Analysis & MAGs
```bash
scripts/BUSCO
```
High-quality bins (MAGs) were selected based on >5% completeness
Quality Control (BUSCO) evaluates the completeness and contamination of the recovered bins 

## 5. Taxonomic Classification of MAGs
Taxonomic classification of MAGs using GTDB-Tk.

```bash
gtdbtk classify_wf \
--genome_dir maxbin2/ \
--extension fa \
--out_dir gtdbtk_output \
--cpus 16
```

## 6. MAG Abundances Profile 
```bash
scripts/cover_M
```
Calculate the real relative abundance of each recovered MAG across all samples. Using CoverM, raw reads are mapped back to the consolidated MAGs database.

## 7. Functional annotation
```bash
scripts/eggNOG_mapper
```
Decodes the metabolic potential of the high-quality MAGs. Using eggNOG-mapper, proteins predicted from your genomes are annotated.

## 8. Genome-Scale Metabolic Modeling
- GEM reconstruction
Runs CarveMe on each MAG nucleotide FASTA to reconstruct genome-scale metabolic models (GEMs) in BiGG namespace using Diamond for functional annotation. Output: one .xml SBML file per genome.  Output: gem_quality_summary.csv and per-genome PDF summaries. Output: one .xml SBML file per genome.
```bash
scripts/01_carveme_reconstruction.sh
```
- GEM quality assessment
Loads each reconstructed GEM with COBRApy, runs FBA, and records reactions, metabolites, genes, growth rate, and solver status per genome. Output: gem_quality_summary.csv and per-genome PDF summaries.
```bash
scripts/02_gem_quality_assessment.py
```
- Community metabolic modeling and exchange network
 Builds a 13-member community with MICOM, runs cooperative tradeoff FBA (fraction=0.5), extracts exchange fluxes, categorizes metabolites (amino acids, carbon sources, cofactors, nucleotides, ions), computes compatibility scores (net contribution 50% + connectivity 30% + diversity 20%), and generates three figures: exchange network for amino acids/carbon (Red 1), exchange network for cofactors/nucleotides/ions (Red 2), and compatibility score barplot. Also generates tabla_interacciones_13genomas.csv and compatibility_scores_13genomas.csv.
```bash
scripts/03_syncom_analysis.py
```
- Inducer/repressor leave-one-out analysis
Performs leave-one-out (LOO) community simulations: removes each genome one at a time, re-runs cooperative tradeoff, and measures the growth rate change of all remaining members. Classifies interactions as inducer (removal decreases growth >5%), repressor (removal increases growth >5%), or neutral. Generates a delta_pct heatmap, an inducer/repressor network, and a summary barplot of how many organisms each genome induces or represses.
```bash
scripts/04_inducer_repressor_analysis.py
```
- Supplementary tables
Assembles a multi-sheet Excel file from all previously generated CSVs: S1 genome summary, S2 all pairwise interactions, S3 producer-consumer pair summary, S4 metabolite summary.
```bash
supplementary/05_supplementary_tables.py
```

# 📊 R Analysis and Data Visualization
The R scripts located in `scripts_R/` were used to analyze the microbial community. The workflow integrates ecological analyses, comparative metagenomics, phylogenetic visualization, and biosynthetic gene cluster characterization to generate publication-quality figures and summary statistics.

The workflow was developed primarily using the following R packages:
```r
install.packages("ggplot2")
install.packages("vegan")
install.packages("tidyverse")
install.packages("patchwork")
install.packages("VennDiagram")
install.packages("UpSetR")
install.packages("ggvenn")
install.packages("RColorBrewer")
install.packages("ape")
```
The analyses included:
### Taxonomic Analysis
Taxonomic profiles were processed to compare microbial community composition across fermented beverages at the genus level. Low-abundance taxa were filtered using a relative abundance threshold of >0.1%, and non-informative assignments such as unclassified or cellular organism labels were excluded to emphasize biologically relevant taxa. 

- Taxonomic Intersections
```bash
R_scripts/upset_venn
```
- Genus-Level Taxonomic Composition
```bash
R_scripts/genus_abundance
R_scripts/family_abundance
```
- Comparative Genus-Level Profile
```bash
R_scripts/comparative_genus
```
### Alpha Diversity
```bash
R_scripts/alpha_diversity
```
Community richness and diversity were summarized using three standard ecological indices: species richness, Shannon diversity, and Simpson diversity. The analysis was performed on genus-level abundance tables after filtering low-quality and non-informative taxonomic assignments.

 ### Beta Diversity
 ```bash
R_scripts/BetaD_NMDS
R_scripts/Umap_BD
```
Beta diversity was assessed using Bray-Curtis dissimilarity on genus-level abundance tables to compare microbial community structure across fermented beverages. Ordination and clustering approaches were used to visualize differences among samples, including UMAP and NMDS-based representations.

### Functional Analysis
 ```bash
R_scripts/BGCs_istribution_heatmap
```
Functional profiles were explored to characterize the metabolic potential of microbial communities and recovered MAGs. Biosynthetic gene clusters (BGCs) identified with antiSMASH were summarized to compare the distribution of specialized metabolite pathways among fermented beverages.

### MAG Characterization 
 ```bash
R_scripts/philogeny_mags
```
Visualizes the phylogenetic placement of metagenome-assembled genomes (MAGs) using GTDB-Tk trees and `ggtree`. The script overlays taxonomic metadata onto the tree, colors tips by genus, and labels MAGs by species to summarize their evolutionary relationships across fermented beverages.

### Comparative and Visualization Analyses
- Integrated Metagenomic Comparison

Integrates MAG taxonomic distribution, BUSCO quality metrics, functional gene abundance, and genus-level intersections into a comparative visualization across fermented beverages.
 ```bash
R_scripts/health_vs_nonhealth
```
- Health vs Non-Health Associated Taxa

Compares the relative abundance of taxa associated with health-related and non-health-related profiles across fermented beverages using log2-transformed abundance values to emphasize low-abundance genera.
 ```bash
R_scripts/BUSCO_TaxDis_Func_G
```

<p align="center">
  <img width="850" alt="Relative Abundance Analysis" src="https://github.com/user-attachments/assets/a91006ae-5685-4abd-8af8-76bac123b1ee" />
</p>
