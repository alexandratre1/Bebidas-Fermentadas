#!/usr/bin/env python3
# =============================================================================
# 02_gem_quality_assessment.py
#
# PURPOSE:
#   Load each reconstructed GEM (SBML .xml), run flux balance analysis (FBA)
#   using COBRApy, and record quality metrics per genome. Also generates a
#   per-genome PDF summary and a global CSV table.
#
# INPUT:
#   CarveMe SBML files at:
#     <BASE>/<genome_id>/<genome_id>_carveme.xml
#
# OUTPUT:
#   - <BASE>/gem_quality_summary.csv       — global metrics table
#   - <BASE>/gem_quality_<genome>.pdf      — per-genome reaction/metabolite plot
#
# DEPENDENCIES:
#   pip install cobra matplotlib pandas
#
# USAGE:
#   conda activate micom_env
#   python3 02_gem_quality_assessment.py
#
# NOTES:
#   - All 13 GEMs should return status "optimal"; any "infeasible" model
#     indicates a reconstruction problem and should be excluded or re-run.
#   - genome1 was flagged as an outlier (only 41 genes) and genome3 showed
#     the highest growth rate (117.0); both were retained but noted in methods.
#   - Exchange reactions are counted as those whose ID starts with "EX_".
#
# AUTHORS: Enrique Pola-Sánchez, Jose Manuel Villalobos-Escobedo
# PROJECT: SynCom design from Mexican fermented beverages
# =============================================================================

import cobra
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
BASE    = "/Users/manuelve/Documents/GEMs/MAGs_new"
GENOMES = [f"genome{i}" for i in list(range(1, 10)) + list(range(10, 14))]

# ------------------------------------------------------------
# Load models, run FBA, collect metrics
# ------------------------------------------------------------
records = []

for genome in GENOMES:
    xml_path = os.path.join(BASE, genome, f"{genome}_carveme.xml")

    if not os.path.isfile(xml_path):
        print(f"[SKIP] {genome}: file not found at {xml_path}")
        continue

    print(f"[LOAD] {genome} ...")
    model  = cobra.io.read_sbml_model(xml_path)
    sol    = model.optimize()

    n_reactions   = len(model.reactions)
    n_metabolites = len(model.metabolites)
    n_genes       = len(model.genes)
    n_exchanges   = sum(1 for r in model.reactions if r.id.startswith("EX_"))
    growth_rate   = sol.objective_value if sol.status == "optimal" else 0.0
    status        = sol.status

    print(f"  reactions={n_reactions}  metabolites={n_metabolites}  "
          f"genes={n_genes}  growth={growth_rate:.4f}  status={status}")

    records.append({
        "genome":       genome,
        "reactions":    n_reactions,
        "metabolites":  n_metabolites,
        "genes":        n_genes,
        "exchange_rxns":n_exchanges,
        "growth_rate":  round(growth_rate, 4),
        "fba_status":   status,
    })

    # ----------------------------------------------------------
    # Per-genome PDF: bar chart of reactions vs metabolites
    # ----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 4))
    categories = ["Reactions", "Metabolites", "Genes", "Exchange rxns"]
    values     = [n_reactions, n_metabolites, n_genes, n_exchanges]
    colors     = ["#2E86C1", "#28B463", "#E67E22", "#8E44AD"]
    bars = ax.bar(categories, values, color=colors, edgecolor="gray", width=0.5)

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                str(val), ha="center", va="bottom", fontsize=10)

    ax.set_title(f"{genome} — GEM quality summary\n"
                 f"Growth rate: {growth_rate:.4f}  |  Status: {status}",
                 fontsize=11, fontweight="bold")
    ax.set_ylabel("Count", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    pdf_path = os.path.join(BASE, f"gem_quality_{genome}.pdf")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {pdf_path}")

# ------------------------------------------------------------
# Save global summary CSV
# ------------------------------------------------------------
df = pd.DataFrame(records)
csv_path = os.path.join(BASE, "gem_quality_summary.csv")
df.to_csv(csv_path, index=False)
print(f"\nSaved global summary: {csv_path}")

# ------------------------------------------------------------
# Print summary table to console
# ------------------------------------------------------------
print("\n=== GEM QUALITY SUMMARY ===")
print(df.to_string(index=False))

# Flag potential outliers
feasible = df[df["fba_status"] == "optimal"]
if len(feasible) < len(df):
    print(f"\n[WARNING] {len(df) - len(feasible)} model(s) did not reach "
          f"optimal status — check those GEMs before community modeling.")
else:
    print(f"\nAll {len(df)} models achieved optimal FBA status.")
