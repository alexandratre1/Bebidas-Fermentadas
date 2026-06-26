#!/usr/bin/env python3
# =============================================================================
# 04_inducer_repressor_analysis.py
#
# PURPOSE:
#   Identify inducer and repressor relationships within the SynCom through
#   leave-one-out (LOO) community simulations. Each genome is removed from
#   the community one at a time, and the resulting change in growth rates
#   of all remaining members is measured relative to the full-community
#   baseline. Interactions are classified as:
#
#     - Inducer  : removing genome X decreases growth of genome Y by >5%
#                  → X promotes Y's growth (positive dependence)
#     - Repressor: removing genome X increases growth of genome Y by >5%
#                  → X suppresses Y's growth (competitive interaction)
#     - Neutral  : |Δ growth| ≤ 5%
#
# INPUT:
#   - CarveMe SBML files at <BASE>/<genome_id>/<genome_id>_carveme.xml
#   - compatibility_scores_13genomas.csv  (from script 03)
#
# OUTPUT:
#   - inducer_repressor_analysis.csv       — full LOO results per pair
#   - inducer_repressor_summary.csv        — summary counts per genome
#   - fig_inducer_repressor_heatmap.pdf    — delta_pct heatmap (all pairs)
#   - fig_inducer_repressor_network.pdf    — directed interaction network
#   - fig_inducer_repressor_summary.pdf    — barplot: induces vs represses
#
# DEPENDENCIES:
#   pip install micom cobra networkx matplotlib pandas numpy
#
# USAGE:
#   conda activate micom_env
#   python3 04_inducer_repressor_analysis.py
#
# KEY PARAMETERS:
#   - tradeoff fraction   = 0.5  (same as community baseline)
#   - effect threshold    = ±5%  (change in growth rate to classify interaction)
#   - compatibility_scores_13genomas.csv used to color nodes
#
# REFERENCES:
#   Zelezniak et al. (2015) PNAS — metabolic dependencies in microbial communities
#   Diener et al. (2020) mSystems — MICOM cooperative tradeoff
#
# AUTHORS: Enrique Pola-Sánchez, Jose Manuel Villalobos-Escobedo
# PROJECT: SynCom design from Mexican fermented beverages
# =============================================================================

import micom
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import FancyArrowPatch
import networkx as nx
from micom import Community

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------
BASE    = "/Users/manuelve/Documents/GEMs/MAGs_new"
GENOMES = [f"genome{i}" for i in list(range(1, 10)) + list(range(10, 14))]

# Effect classification thresholds (% change in growth rate)
EFFECT_THRESHOLD = 5.0   # |Δ%| > 5 → inducer or repressor

# Colors for interaction types
COLOR_EFFECT = {
    "inducer":   "#2ECC71",   # green — removal decreases growth (X supports Y)
    "repressor": "#E74C3C",   # red   — removal increases growth (X competes with Y)
}

# ============================================================
# STEP 1: Run full-community baseline simulation
# ============================================================
manifest_full = pd.DataFrame({
    "id":        GENOMES,
    "file":      [f"{BASE}/{g}/{g}_carveme.xml" for g in GENOMES],
    "abundance": [1 / len(GENOMES)] * len(GENOMES)
})

print("Running full-community baseline...")
com_full     = Community(manifest_full, progress=False)
result_full  = com_full.cooperative_tradeoff(fraction=0.5, fluxes=False, pfba=False)
baseline     = result_full.members["growth_rate"].to_dict()

print("Baseline growth rates:")
for g, r in baseline.items():
    print(f"  {g}: {r:.4f}")

# ============================================================
# STEP 2: Leave-one-out simulations
# ============================================================
print("\nRunning leave-one-out simulations...")
loo_results = {}

for removed in GENOMES:
    remaining = [g for g in GENOMES if g != removed]
    n         = len(remaining)
    manifest_loo = pd.DataFrame({
        "id":        remaining,
        "file":      [f"{BASE}/{g}/{g}_carveme.xml" for g in remaining],
        "abundance": [1 / n] * n
    })
    print(f"  Removing {removed} ({n} remaining)...")
    try:
        com_loo    = Community(manifest_loo, progress=False)
        result_loo = com_loo.cooperative_tradeoff(fraction=0.5, fluxes=False, pfba=False)
        loo_results[removed] = result_loo.members["growth_rate"].to_dict()
    except Exception as e:
        print(f"    ERROR: {e}")
        loo_results[removed] = {g: np.nan for g in remaining}

# ============================================================
# STEP 3: Classify interactions (inducer / repressor / neutral)
# ============================================================
rows = []
for removed, growth_dict in loo_results.items():
    for target, new_rate in growth_dict.items():
        baseline_rate = baseline.get(target, np.nan)
        if np.isnan(baseline_rate) or baseline_rate == 0:
            delta     = np.nan
            delta_pct = np.nan
        else:
            delta     = new_rate - baseline_rate
            delta_pct = (delta / baseline_rate) * 100

        # Classify based on direction and magnitude of growth rate change
        if np.isnan(delta_pct):
            effect = "unknown"
        elif delta_pct < -EFFECT_THRESHOLD:
            # Removing 'removed' decreased 'target' growth → 'removed' is an inducer of 'target'
            effect = "inducer"
        elif delta_pct > EFFECT_THRESHOLD:
            # Removing 'removed' increased 'target' growth → 'removed' is a repressor of 'target'
            effect = "repressor"
        else:
            effect = "neutral"

        rows.append({
            "removed":       removed,
            "target":        target,
            "baseline_rate": round(baseline_rate, 6),
            "new_rate":      round(new_rate, 6),
            "delta":         round(delta, 6) if not np.isnan(delta) else np.nan,
            "delta_pct":     round(delta_pct, 4) if not np.isnan(delta_pct) else np.nan,
            "effect":        effect,
        })

df = pd.DataFrame(rows)
df.to_csv(f"{BASE}/inducer_repressor_analysis.csv", index=False)
print(f"\nSaved: inducer_repressor_analysis.csv  ({len(df)} rows)")

# Print significant interactions
print("\n=== INDUCTORS AND REPRESSORS ===")
df_sig = df[df["effect"] != "neutral"].sort_values("delta_pct")
print(df_sig[["removed", "target", "baseline_rate", "new_rate",
              "delta_pct", "effect"]].to_string(index=False))

# ============================================================
# FIGURE 1: Heatmap of Δ growth rate (%) for all LOO pairs
# ============================================================
pivot = df.pivot(index="target", columns="removed", values="delta_pct")

fig1, ax1 = plt.subplots(figsize=(14, 10))
vals = pivot.values[~np.isnan(pivot.values)]
vmax = max(abs(vals.max()), abs(vals.min()))

im = ax1.imshow(pivot.values, cmap="RdBu", vmin=-vmax, vmax=vmax, aspect="auto")
ax1.set_xticks(range(len(pivot.columns)))
ax1.set_yticks(range(len(pivot.index)))
ax1.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=10)
ax1.set_yticklabels(pivot.index, fontsize=10)
ax1.set_xlabel("Removed genome", fontsize=12)
ax1.set_ylabel("Target genome", fontsize=12)
ax1.set_title(
    "Effect of Genome Removal on Growth Rates\n(% change from baseline)",
    fontsize=13, fontweight="bold")

# Annotate each cell with the numeric value
for i in range(len(pivot.index)):
    for j in range(len(pivot.columns)):
        val = pivot.values[i, j]
        if not np.isnan(val):
            ax1.text(j, i, f"{val:.1f}", ha="center", va="center",
                     fontsize=7,
                     color="black" if abs(val) < vmax * 0.6 else "white")

plt.colorbar(im, ax=ax1, label="Δ Growth Rate (%)")
plt.tight_layout()
plt.savefig(f"{BASE}/fig_inducer_repressor_heatmap.pdf", bbox_inches="tight")
plt.close()
print("Saved: fig_inducer_repressor_heatmap.pdf")

# ============================================================
# FIGURE 2: Directed inducer/repressor network
# ============================================================
# Load compatibility scores for node coloring
df_scores   = pd.read_csv(f"{BASE}/compatibility_scores_13genomas.csv", index_col=0)
norm_nodes  = Normalize(vmin=df_scores["compatibility_score"].min(),
                        vmax=df_scores["compatibility_score"].max())
cmap_nodes  = plt.cm.RdBu_r
node_colors = [cmap_nodes(norm_nodes(df_scores.loc[n, "compatibility_score"]))
               for n in GENOMES]

# Build directed graph (one edge per significant inducer/repressor pair)
G = nx.DiGraph()
G.add_nodes_from(GENOMES)

df_sig_net = df[df["effect"].isin(["inducer", "repressor"])]
for _, row in df_sig_net.iterrows():
    u, v   = row["removed"], row["target"]
    effect = row["effect"]
    w      = abs(row["delta_pct"])
    if G.has_edge(u, v):
        # Keep the strongest interaction if multiple exist
        if w > abs(G[u][v]["delta_pct"]):
            G[u][v].update({"weight": w, "effect": effect,
                            "delta_pct": row["delta_pct"]})
    else:
        G.add_edge(u, v, weight=w, effect=effect, delta_pct=row["delta_pct"])

weights = [G[u][v]["weight"] for u, v in G.edges()] if G.edges() else [1]
max_w   = max(weights) if weights else 1

fig2, ax2 = plt.subplots(figsize=(14, 12))
pos = nx.circular_layout(G)
ax2.set_xlim(-1.5, 1.5)
ax2.set_ylim(-1.5, 1.5)

nx.draw_networkx_nodes(G, pos, node_size=300,
                       node_color=node_colors, edgecolors="#333333",
                       linewidths=2, ax=ax2)
nx.draw_networkx_labels(G, pos, {n: n for n in G.nodes()},
                        font_size=9, font_weight="bold", ax=ax2)

for u, v in G.edges():
    w      = G[u][v]["weight"]
    lw     = max(1.5, 5 * w / max_w)
    ec     = COLOR_EFFECT[G[u][v]["effect"]]
    x1, y1 = pos[u]
    x2, y2 = pos[v]
    arrow  = FancyArrowPatch(
        posA=(x1, y1), posB=(x2, y2),
        arrowstyle="-|>", color=ec, linewidth=lw,
        mutation_scale=18, connectionstyle="arc3,rad=0.2",
        shrinkA=12, shrinkB=12, zorder=3)
    ax2.add_patch(arrow)

    # Label with the delta_pct value
    mx = (x1 + x2) / 2
    my = (y1 + y2) / 2
    dx, dy = x2 - x1, y2 - y1
    nv = (dx**2 + dy**2)**0.5 + 1e-6
    ox, oy = -dy/nv * 0.1, dx/nv * 0.1
    ax2.text(mx + ox, my + oy, f"{G[u][v]['delta_pct']:.1f}%",
             fontsize=6, ha="center", va="center",
             bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=ec,
                       alpha=0.9, linewidth=0.8), zorder=5)

patches = [
    mpatches.Patch(color="#2ECC71", label="Inducer (removal decreases growth)"),
    mpatches.Patch(color="#E74C3C", label="Repressor (removal increases growth)"),
]
ax2.legend(handles=patches, loc="lower left", fontsize=9, framealpha=0.9)

sm = ScalarMappable(cmap=cmap_nodes, norm=norm_nodes)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax2, shrink=0.4, pad=0.02)
cbar.set_label("Compatibility Score", fontsize=9)

ax2.set_title("Inducer / Repressor Network\n(green = inducer, red = repressor)",
              fontsize=13, fontweight="bold")
ax2.axis("off")
plt.tight_layout()
plt.savefig(f"{BASE}/fig_inducer_repressor_network.pdf", bbox_inches="tight")
plt.close()
print("Saved: fig_inducer_repressor_network.pdf")

# ============================================================
# FIGURE 3: Summary barplot — how many organisms each genome
#           induces/represses and is induced/repressed by
# ============================================================
summary_rows = []
for g in GENOMES:
    summary_rows.append({
        "genome":       g,
        "induces":      len(df[(df["removed"] == g) & (df["effect"] == "inducer")]),
        "represses":    len(df[(df["removed"] == g) & (df["effect"] == "repressor")]),
        "is_induced":   len(df[(df["target"] == g) & (df["effect"] == "inducer")]),
        "is_repressed": len(df[(df["target"] == g) & (df["effect"] == "repressor")]),
    })

df_sum = pd.DataFrame(summary_rows)
df_sum.to_csv(f"{BASE}/inducer_repressor_summary.csv", index=False)
print("Saved: inducer_repressor_summary.csv")

fig3, axes = plt.subplots(1, 2, figsize=(16, 8))
x = np.arange(len(GENOMES))
w = 0.35

# Left panel: how many organisms each genome induces / represses
ax3 = axes[0]
ax3.bar(x - w/2, df_sum["induces"],   width=w, color="#2ECC71",
        edgecolor="gray", label="Induces others")
ax3.bar(x + w/2, df_sum["represses"], width=w, color="#E74C3C",
        edgecolor="gray", label="Represses others")
ax3.set_xticks(x)
ax3.set_xticklabels(GENOMES, rotation=45, ha="right", fontsize=9)
ax3.set_ylabel("Number of genomes affected", fontsize=11)
ax3.set_title("How many genomes each\norganism induces or represses",
              fontsize=12, fontweight="bold")
ax3.legend(fontsize=9)
ax3.spines["top"].set_visible(False)
ax3.spines["right"].set_visible(False)

# Right panel: how many organisms induce / repress each genome
ax4 = axes[1]
ax4.bar(x - w/2, df_sum["is_induced"],   width=w, color="#2ECC71",
        edgecolor="gray", label="Induced by others")
ax4.bar(x + w/2, df_sum["is_repressed"], width=w, color="#E74C3C",
        edgecolor="gray", label="Repressed by others")
ax4.set_xticks(x)
ax4.set_xticklabels(GENOMES, rotation=45, ha="right", fontsize=9)
ax4.set_ylabel("Number of genomes affecting it", fontsize=11)
ax4.set_title("How many genomes\ninduce or repress each organism",
              fontsize=12, fontweight="bold")
ax4.legend(fontsize=9)
ax4.spines["top"].set_visible(False)
ax4.spines["right"].set_visible(False)

plt.suptitle("Inducer / Repressor Summary — 13-Genome SynCom",
             fontsize=14, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig(f"{BASE}/fig_inducer_repressor_summary.pdf", bbox_inches="tight")
plt.close()
print("Saved: fig_inducer_repressor_summary.pdf")
