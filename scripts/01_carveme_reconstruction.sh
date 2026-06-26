#!/usr/bin/env bash
# =============================================================================
# 01_carveme_reconstruction.sh
#
# PURPOSE:
#   Reconstruct genome-scale metabolic models (GEMs) from metagenome-assembled
#   genome (MAG) nucleotide FASTA files using CarveMe v1.6.6.
#
#   CarveMe internally runs:
#     1. Prodigal  — gene prediction from nucleotide sequence
#     2. Diamond   — functional annotation against the BiGG reference database
#     3. Model carving + gap-filling → SBML output in BiGG namespace
#
# INPUT:
#   Nucleotide FASTA files (.fasta) located at:
#     <BASE>/<genome_id>/<genome_id>.fasta
#
# OUTPUT:
#   One SBML file (.xml) per genome at:
#     <BASE>/<genome_id>/<genome_id>_carveme.xml
#
# DEPENDENCIES:
#   - CarveMe v1.6.6  (pip install carveme)
#   - Diamond v2.1.8  (conda install -c bioconda diamond)
#   - Miniforge / conda environment: micom_env (Python 3.11)
#
# USAGE:
#   conda activate micom_env
#   bash 01_carveme_reconstruction.sh
#
# NOTES:
#   - The --dna flag is required when the input is a nucleotide FASTA
#     (as opposed to a protein FASTA). CarveMe triggers Prodigal internally.
#   - All 13 genomes produced feasible GEMs (solver status: optimal).
#   - Estimated runtime: ~2–5 min per genome on a standard laptop.
#
# AUTHORS: Enrique Pola-Sánchez, Jose Manuel Villalobos-Escobedo
# PROJECT: SynCom design from Mexican fermented beverages
# =============================================================================

set -euo pipefail

# ------------------------------------------------------------
# Configuration — adjust BASE path to match your directory
# ------------------------------------------------------------
BASE="/Users/manuelve/Documents/GEMs/MAGs_new"

GENOMES=(
    genome1  genome2  genome3  genome4  genome5
    genome6  genome7  genome8  genome9  genome10
    genome11 genome12 genome13
)

# ------------------------------------------------------------
# Reconstruct GEMs — one per genome
# ------------------------------------------------------------
echo "Starting CarveMe GEM reconstruction for ${#GENOMES[@]} genomes..."
echo "Base directory: ${BASE}"
echo ""

for GENOME in "${GENOMES[@]}"; do
    FASTA="${BASE}/${GENOME}/${GENOME}.fasta"
    OUTPUT="${BASE}/${GENOME}/${GENOME}_carveme.xml"

    if [[ ! -f "${FASTA}" ]]; then
        echo "  [SKIP] ${GENOME}: FASTA not found at ${FASTA}"
        continue
    fi

    if [[ -f "${OUTPUT}" ]]; then
        echo "  [SKIP] ${GENOME}: output already exists at ${OUTPUT}"
        continue
    fi

    echo "  [RUN]  ${GENOME} ..."
    carve "${FASTA}" \
          --dna \
          -o "${OUTPUT}"

    echo "  [DONE] ${GENOME} → ${OUTPUT}"
done

echo ""
echo "All reconstructions complete."
echo "Output files:"
for GENOME in "${GENOMES[@]}"; do
    XML="${BASE}/${GENOME}/${GENOME}_carveme.xml"
    if [[ -f "${XML}" ]]; then
        echo "  OK  ${XML}"
    else
        echo "  MISSING  ${XML}"
    fi
done
