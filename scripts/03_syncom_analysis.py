#!/usr/bin/env python3
# =============================================================================
# 03_syncom_analysis.py
#
# PURPOSE:
#   Build a 13-member synthetic community (SynCom) using MICOM, run
#   cooperative tradeoff FBA, extract metabolic exchange fluxes, compute
#   per-genome compatibility scores, and generate metabolic exchange network
#   visualizations.
#
# ANALYSIS STEPS:
#   1. Assemble community manifest (equal abundances, 1/n per member)
#   2. Run cooperative tradeoff FBA (fraction=0.5, fluxes=True)
#   3. Extract and filter exchange fluxes (exclude stoichiometric artifacts)
#   4. Categorize exchanged metabolites into functional groups
#   5. Compute compatibility scores per genome:
#        net_contribution (50%) + connectivity (30%) + diversity (20%)
#   6. Plot Red 1: amino acid & carbon source exchange network
#   7. Plot Red 2: cofactor / nucleotide / ion exchange network
#   8. Plot Figure 3: compatibility score barplot per genome
#
# INPUT:
#   CarveMe SBML files at <BASE>/<genome_id>/<genome_id>_carveme.xml
#
# OUTPUT:
#   - tabla_interacciones_13genomas.csv   — all pairwise flux interactions
#   - compatibility_scores_13genomas.csv  — per-genome compatibility metrics
#   - red1_aminoacidos_carbono.pdf        — amino acid & carbon exchange network
#   - red2_cofactores_nucleotidos_iones.pdf — cofactor/nucleotide/ion network
#   - fig3_compatibility_scores.pdf       — compatibility score barplot
#
# DEPENDENCIES:
#   pip install micom cobra networkx matplotlib pandas
#
# USAGE:
#   conda activate micom_env
#   python3 03_syncom_analysis.py
#
# KEY PARAMETERS:
#   - tradeoff fraction = 0.5  (balances individual vs community growth)
#   - significance_pct  > 0.1  (minimum % of total flux to include interaction)
#   - compatibility threshold = 0.35  (minimum score to label genome as compatible)
#   - Excluded metabolites: h2o, h, co2, pi, nh4, o2  (stoichiometric artifacts)
#
# AUTHORS: Enrique Pola-Sánchez, Jose Manuel Villalobos-Escobedo
# PROJECT: SynCom design from Mexican fermented beverages
# =============================================================================

import micom
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import FancyArrowPatch
from micom import Community

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
BASE    = "/Users/manuelve/Documents/GEMs/MAGs_new"
GENOMES = [f"genome{i}" for i in list(range(1, 10)) + list(range(10, 14))]

# Metabolites excluded from exchange analysis (stoichiometric artifacts /
# inorganic compounds that do not represent meaningful cross-feeding)
EXCLUDE = {"h2o", "h", "co2", "pi", "nh4", "o2"}

# Significance threshold: minimum % of a genome's total flux for an
# exchange to be considered biologically relevant
SIGNIFICANCE_THRESHOLD = 0.1

# Compatibility score threshold: genomes scoring >= this value are
# considered "compatible" SynCom members
COMPATIBILITY_THRESHOLD = 0.35

# ------------------------------------------------------------
# Metabolite display names (BiGG ID → human-readable)
# ------------------------------------------------------------
MET_FULL_NAMES = {
    "ac": "Acetate", "pep": "Phosphoenolpyruvate",
    "isobuta": "Isobutyrate", "val__L": "L-Valine",
    "hom__L": "L-Homoserine", "orn": "L-Ornithine",
    "dhap": "Dihydroxyacetone phosphate", "g6p": "Glucose-6-phosphate",
    "g6p_B": "Glucose-6-phosphate (B)", "f6p": "Fructose-6-phosphate",
    "thm": "Thiamine", "coa": "Coenzyme A", "etoh": "Ethanol",
    "succ": "Succinate", "pyr": "Pyruvate", "glu__L": "L-Glutamate",
    "gln__L": "L-Glutamine", "asp__L": "L-Aspartate",
    "leu__L": "L-Leucine", "ile__L": "L-Isoleucine",
    "lys__L": "L-Lysine", "arg__L": "L-Arginine",
    "ser__L": "L-Serine", "thr__L": "L-Threonine",
    "met__L": "L-Methionine", "pro__L": "L-Proline",
    "asn__L": "L-Asparagine", "LalaDgluMdap": "L-Ala-D-Glu-mDAP",
    "mal__L": "L-Malate", "lac__D": "D-Lactate", "lac__L": "L-Lactate",
    "amp": "AMP", "gmp": "GMP", "cmp": "CMP", "ump": "UMP",
    "adn": "Adenosine", "gsn": "Guanosine",
    "fe2": "Iron (Fe2+)", "fe3": "Iron (Fe3+)",
    "mg2": "Magnesium (Mg2+)", "so4": "Sulfate",
    "glc__D": "D-Glucose", "fru": "Fructose", "man": "Mannose",
    "gal": "Galactose", "tre": "Trehalose", "xyl__D": "D-Xylose",
    "for": "Formate", "but": "Butyrate", "prop": "Propionate",
    "glyc": "Glycerol", "gly": "Glycine", "trp__L": "L-Tryptophan",
    "phe__L": "L-Phenylalanine", "tyr__L": "L-Tyrosine",
    "his__L": "L-Histidine", "cys__L": "L-Cysteine",
    "ala__L": "L-Alanine", "fol": "Folate", "ribflv": "Riboflavin",
    "nac": "Nicotinate", "pnto__R": "Pantothenate",
    "btn": "Biotin", "k": "Potassium (K+)", "na1": "Sodium (Na+)",
    "cl": "Chloride (Cl-)", "ca2": "Calcium (Ca2+)",
    "zn2": "Zinc (Zn2+)", "mn2": "Manganese (Mn2+)",
    "cbl1": "Cob(I)alamin (B12)", "adocbl": "Adenosylcobalamin (B12)"
}

# Metabolite functional categories (BiGG IDs mapped to category labels)
CATEGORIES = {
    "Amino acids":    ["val__L", "leu__L", "ile__L", "glu__L", "gln__L",
                       "asp__L", "asn__L", "orn", "hom__L", "ser__L",
                       "pro__L", "arg__L", "lys__L", "thr__L", "met__L",
                       "gly", "trp__L", "phe__L", "tyr__L", "his__L",
                       "cys__L", "ala__L"],
    "Carbon sources": ["ac", "succ", "mal__L", "pyr", "etoh", "lac__D",
                       "lac__L", "isobuta", "pep", "f6p", "g6p", "dhap",
                       "for", "but", "prop", "glyc", "glc__D", "fru",
                       "man", "gal", "tre", "xyl__D"],
    "Cofactors":      ["fol", "ribflv", "nac", "pnto__R", "btn", "thm",
                       "fe2", "fe3", "mg2", "so4", "zn2", "mn2",
                       "cbl1", "adocbl"],
    "Nucleotides":    ["amp", "gmp", "cmp", "ump", "adn", "gsn"],
    "Ions":           ["k", "na1", "cl", "ca2"]
}

# Colors for each network visualization
COLOR_NET1 = {"Amino acids": "#E74C3C", "Carbon sources": "#2ECC71"}
COLOR_NET2 = {"Cofactors": "#3498DB", "Nucleotides": "#9B59B6", "Ions": "#F39C12"}

# Colormap for node compatibility scores
CMAP_NODES = plt.cm.RdBu_r


def get_full_name(met_id):
    """Return human-readable metabolite name from BiGG ID."""
    return MET_FULL_NAMES.get(met_id, met_id)


def get_category(met_id):
    """Assign a metabolite to its functional category based on BiGG ID."""
    for cat, mets in CATEGORIES.items():
        if any(m in met_id for m in mets):
            return cat
    return "Other"


# ============================================================
# STEP 1: Build community and run cooperative tradeoff FBA
# ============================================================
manifest = pd.DataFrame({
    "id":        GENOMES,
    "file":      [f"{BASE}/{g}/{g}_carveme.xml" for g in GENOMES],
    "abundance": [1 / len(GENOMES)] * len(GENOMES)
})

print("Building community model...")
com = Community(manifest, progress=True)

print("Running cooperative tradeoff FBA (fraction=0.5)...")
result = com.cooperative_tradeoff(fraction=0.5, fluxes=True, pfba=False)
fluxes = result.fluxes

# Total absolute flux per genome (used to compute significance_pct)
total_flux = fluxes.abs().sum(axis=1)

# ============================================================
# STEP 2: Extract pairwise exchange interactions
# ============================================================
rows = []
for met in fluxes.columns:
    if "EX_" not in met:
        continue
    met_id   = met.replace("EX_", "").replace("_e", "")
    if met_id in EXCLUDE:
        continue
    met_name = get_full_name(met_id)
    category = get_category(met_id)

    # Identify producers (positive flux = secretion) and
    # consumers (negative flux = uptake) for this metabolite
    producers, consumers = [], []
    for genome in GENOMES:
        flux = fluxes.loc[genome, met]
        if flux > 0.0001:
            producers.append((genome, flux))
        elif flux < -0.0001:
            consumers.append((genome, abs(flux)))

    # Record each producer→consumer pair as a directed interaction
    for prod, f_prod in producers:
        for cons, f_cons in consumers:
            if prod != cons:
                rows.append({
                    "producer":         prod,
                    "consumer":         cons,
                    "metabolite_id":    met_id,
                    "metabolite_name":  met_name,
                    "category":         category,
                    "flux_producer":    round(f_prod, 6),
                    "flux_consumer":    round(f_cons, 6),
                    "exchange_flux":    round((f_prod + f_cons) / 2, 6),
                    # significance_pct: fraction of the producer's total flux
                    # dedicated to exporting this metabolite (not a p-value)
                    "significance_pct": round(f_prod / total_flux[prod] * 100, 4),
                })

df_all = pd.DataFrame(rows).sort_values("significance_pct", ascending=False)
df_all.to_csv(f"{BASE}/tabla_interacciones_13genomas.csv", index=False)
print(f"Saved: tabla_interacciones_13genomas.csv  ({len(df_all)} interactions)")

# Apply significance filter
df = df_all[df_all["significance_pct"] > SIGNIFICANCE_THRESHOLD]
print(f"Interactions after significance filter (>{SIGNIFICANCE_THRESHOLD}%): {len(df)}")

# ============================================================
# STEP 3: Compute per-genome compatibility scores
# ============================================================
scores = {}
for genome in GENOMES:
    exported       = df[df["producer"] == genome]["exchange_flux"].sum()
    imported       = df[df["consumer"] == genome]["exchange_flux"].sum()
    partners_out   = df[df["producer"] == genome]["consumer"].nunique()
    partners_in    = df[df["consumer"] == genome]["producer"].nunique()
    mets_exported  = df[df["producer"] == genome]["metabolite_id"].nunique()
    mets_imported  = df[df["consumer"] == genome]["metabolite_id"].nunique()

    scores[genome] = {
        "exported_flux":    round(exported, 4),
        "imported_flux":    round(imported, 4),
        # Net contribution: positive = net donor, negative = net receiver
        "net_contribution": round(exported - imported, 4),
        # Connectivity: total number of unique exchange partners
        "connectivity":     partners_out + partners_in,
        # Diversity: total unique metabolite types exchanged
        "diversity":        mets_exported + mets_imported,
        "growth_rate":      round(result.members.loc[genome, "growth_rate"], 6),
    }

df_scores = pd.DataFrame(scores).T

# Min-max normalize each component before weighted combination
for col in ["net_contribution", "connectivity", "diversity"]:
    min_v = df_scores[col].min()
    max_v = df_scores[col].max()
    df_scores[f"{col}_norm"] = (
        (df_scores[col] - min_v) / (max_v - min_v)
        if max_v > min_v else 0.5
    )

# Weighted compatibility score (weights reflect biological priority)
df_scores["compatibility_score"] = (
    df_scores["net_contribution_norm"] * 0.50 +   # main donor capacity
    df_scores["connectivity_norm"]     * 0.30 +   # integration in network
    df_scores["diversity_norm"]        * 0.20     # metabolic versatility
)
df_scores["compatible"] = df_scores["compatibility_score"] >= COMPATIBILITY_THRESHOLD
df_scores.to_csv(f"{BASE}/compatibility_scores_13genomas.csv")
print(f"Saved: compatibility_scores_13genomas.csv")

print("\n=== COMPATIBILITY SCORES ===")
print(df_scores[["net_contribution", "connectivity", "diversity",
                  "compatibility_score", "compatible"]].to_string())

# Node colors for network figures (mapped from compatibility score)
norm_nodes  = Normalize(vmin=df_scores["compatibility_score"].min(),
                        vmax=df_scores["compatibility_score"].max())
node_colors = {n: CMAP_NODES(norm_nodes(df_scores.loc[n, "compatibility_score"]))
               for n in GENOMES}


# ============================================================
# Helper: build directed NetworkX graph from filtered df
# ============================================================
def build_graph(df_filtered):
    """
    Construct a directed graph where each edge (producer → consumer)
    accumulates total exchange flux and records all exchanged metabolites.
    """
    G = nx.DiGraph()
    G.add_nodes_from(GENOMES)
    for _, row in df_filtered.iterrows():
        prod, cons = row["producer"], row["consumer"]
        met_name   = row["metabolite_name"]
        cat        = row["category"]
        f_prod     = row["flux_producer"]
        if G.has_edge(prod, cons):
            G[prod][cons]["weight"]     += f_prod
            G[prod][cons]["metabolites"].append(met_name)
            G[prod][cons]["categories"].add(cat)
        else:
            G.add_edge(prod, cons, weight=f_prod,
                       metabolites=[met_name], categories={cat})
    return G


# ============================================================
# Helper: draw metabolic exchange network
# ============================================================
def draw_network(ax, G, title, node_colors, color_cat):
    """
    Draw a circular exchange network with:
      - Node color = compatibility score (RdBu colormap)
      - Edge color = metabolite category
      - Edge width = proportional to exchange flux
      - Label = top exchanged metabolite per edge
    """
    pos = nx.circular_layout(G)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)

    colors = [node_colors[n] for n in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_size=300,
                           node_color=colors, edgecolors="#333333",
                           linewidths=2, ax=ax)
    nx.draw_networkx_labels(G, pos, {n: n for n in G.nodes()},
                            font_size=9, font_weight="bold", ax=ax)

    weights = [G[u][v]["weight"] for u, v in G.edges()] if G.edges() else [1]
    max_w   = max(weights) if weights else 1

    for u, v in G.edges():
        w        = G[u][v]["weight"]
        lw       = max(1.5, 5 * w / max_w)
        x1, y1   = pos[u]
        x2, y2   = pos[v]
        dx, dy   = x2 - x1, y2 - y1
        nv       = (dx**2 + dy**2)**0.5 + 1e-6
        cat_list = sorted([c for c in G[u][v]["categories"] if c in color_cat])

        if not cat_list:
            continue

        # Draw one arrow per category with increasing arc radius
        for i, cat in enumerate(cat_list):
            ec  = color_cat[cat]
            rad = 0.15 + i * 0.15
            arrow = FancyArrowPatch(
                posA=(x1, y1), posB=(x2, y2),
                arrowstyle="-|>", color=ec, linewidth=lw,
                mutation_scale=18,
                connectionstyle=f"arc3,rad={rad}",
                shrinkA=12, shrinkB=12, zorder=2 + i
            )
            ax.add_patch(arrow)

        # Label with the top metabolite for this edge
        main_cat = cat_list[0]
        label    = G[u][v]["metabolites"][0]
        mx = (x1 + x2) / 2
        my = (y1 + y2) / 2
        ox, oy = -dy/nv * 0.13, dx/nv * 0.13
        ax.text(mx + ox, my + oy, label, fontsize=6,
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white",
                          ec=color_cat.get(main_cat, "#95A5A6"),
                          alpha=0.9, linewidth=0.8), zorder=10)

    # Legend: metabolite categories
    patches = [mpatches.Patch(color=c, label=l) for l, c in color_cat.items()]
    ax.legend(handles=patches, loc="lower left", fontsize=8,
              title="Metabolite category", title_fontsize=9, framealpha=0.9)

    # Colorbar: compatibility score
    sm = ScalarMappable(cmap=CMAP_NODES, norm=norm_nodes)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.4, pad=0.02)
    cbar.set_label("Compatibility Score", fontsize=9)
    cbar.ax.axhline(y=norm_nodes(COMPATIBILITY_THRESHOLD),
                    color="black", linewidth=1.5, linestyle="--")

    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.axis("off")


# ============================================================
# FIGURE 1 (Red 1): Amino acid & carbon source exchanges
# ============================================================
df_red1 = df[df["category"].isin({"Amino acids", "Carbon sources"})]
G1 = build_graph(df_red1)
print(f"\nNetwork 1 — amino acids + carbon: {len(df_red1)} rows, "
      f"{G1.number_of_edges()} edges")

fig1, ax1 = plt.subplots(figsize=(14, 12))
draw_network(ax1, G1,
             "Amino Acid & Carbon Source Exchanges\n"
             f"(significance > {SIGNIFICANCE_THRESHOLD}%)",
             node_colors, COLOR_NET1)
plt.tight_layout()
plt.savefig(f"{BASE}/red1_aminoacidos_carbono.pdf", bbox_inches="tight")
plt.close()
print("Saved: red1_aminoacidos_carbono.pdf")

# ============================================================
# FIGURE 2 (Red 2): Cofactor / nucleotide / ion exchanges
# ============================================================
df_red2 = df[df["category"].isin({"Cofactors", "Nucleotides", "Ions"})]
G2 = build_graph(df_red2)
print(f"\nNetwork 2 — cofactors/nucleotides/ions: {len(df_red2)} rows, "
      f"{G2.number_of_edges()} edges")

fig2, ax2 = plt.subplots(figsize=(14, 12))
draw_network(ax2, G2,
             "Cofactor / Nucleotide / Ion Exchanges\n"
             f"(significance > {SIGNIFICANCE_THRESHOLD}%)",
             node_colors, COLOR_NET2)
plt.tight_layout()
plt.savefig(f"{BASE}/red2_cofactores_nucleotidos_iones.pdf", bbox_inches="tight")
plt.close()
print("Saved: red2_cofactores_nucleotidos_iones.pdf")

# ============================================================
# FIGURE 3: Compatibility score barplot
# ============================================================
fig3, ax3 = plt.subplots(figsize=(10, 8))
scores_vals = df_scores.loc[GENOMES, "compatibility_score"]
colors_bar  = [CMAP_NODES(norm_nodes(s)) for s in scores_vals]
bars = ax3.barh(GENOMES, scores_vals, color=colors_bar, edgecolor="gray")
ax3.axvline(x=COMPATIBILITY_THRESHOLD, color="black", linewidth=1.5,
            linestyle="--", label=f"Threshold ({COMPATIBILITY_THRESHOLD})")
ax3.set_xlabel("Compatibility Score (0–1)", fontsize=12)
ax3.set_title(
    "Compatibility Score per Genome\n"
    "(net contribution 50% + connectivity 30% + diversity 20%)",
    fontsize=12, fontweight="bold")
ax3.legend(fontsize=10)

for bar, genome in zip(bars, GENOMES):
    score = df_scores.loc[genome, "compatibility_score"]
    ax3.text(score + 0.01, bar.get_y() + bar.get_height() / 2,
             f"{score:.3f}", va="center", fontsize=10)

ax3.set_xlim(0, 1.15)
ax3.invert_yaxis()
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)

sm3 = ScalarMappable(cmap=CMAP_NODES, norm=norm_nodes)
sm3.set_array([])
cbar3 = plt.colorbar(sm3, ax=ax3, shrink=0.5, pad=0.02)
cbar3.set_label("Compatibility Score", fontsize=10)
cbar3.ax.axhline(y=norm_nodes(COMPATIBILITY_THRESHOLD),
                 color="black", linewidth=1.5, linestyle="--")

plt.tight_layout()
plt.savefig(f"{BASE}/fig3_compatibility_scores.pdf", bbox_inches="tight")
plt.close()
print("Saved: fig3_compatibility_scores.pdf")
