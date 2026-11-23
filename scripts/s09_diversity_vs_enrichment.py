#!/usr/bin/env python3
"""
R1-Major8: Mash-distance matrix + relation between within-group diversity and the number of enriched k-mers (Fig. S[X], Table S[X]).
Inputs (per molecule project):  output/diversity_mash_matrix_k<k>.csv and output/diversity_mash_k<k>.csv  (from s03_diversity_kmer_mash.py)
                                output/enriched_<k>mers_stats.csv   (MAIN full-data run; same file as used for Fig. 1)
usage:  python scripts/s09_diversity_vs_enrichment.py -k 15 --dna project_rebuttal --rna project_rebuttal_RNA [--out FigS_diversity]
Output: <out>.svg (top: Mash matrices; bottom: within-group diversity vs enriched k-mers, Spearman rho with EXACT permutation p), <out>_table.csv
"""
import argparse, itertools, numpy as np, pandas as pd
from scipy.stats import rankdata
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ALIAS = {"R": "Rare"}; ORDER = ["A", "B", "C", "D", "Rare", "CRF01.AE", "CRF02.AG"]
ap = argparse.ArgumentParser(); ap.add_argument("-k", type=int, default=15); ap.add_argument("--dna", required=True); ap.add_argument("--rna", required=True); ap.add_argument("--out", default="FigS_diversity")
a = ap.parse_args()
def exact_spearman(x, y):
    rx, ry = rankdata(x), rankdata(y); rho = np.corrcoef(rx, ry)[0, 1]
    perms = [np.corrcoef(rx, ry[list(p)])[0, 1] for p in itertools.permutations(range(len(y)))]      # all n! permutations (n = 7: 5,040)
    return rho, np.mean(np.abs(perms) >= abs(rho) - 1e-12)
def load(proj):
    d = pd.read_csv(f"{proj}/output/diversity_mash_k{a.k}.csv", index_col=0); d.index = [ALIAS.get(i, i) for i in d.index]
    s = pd.read_csv(f"{proj}/output/enriched_{a.k}mers_stats.csv"); s = s[s.group.str.endswith("_enriched")]
    n = s.group.str.replace("_enriched", "").map(lambda g: ALIAS.get(g, g)).value_counts()
    d["enriched_kmers"] = n.reindex(d.index); d = d.reindex([g for g in ORDER if g in d.index])
    try: M = pd.read_csv(f"{proj}/output/diversity_mash_matrix_k{a.k}.csv", index_col=0); M.index = [ALIAS.get(i, i) for i in M.index]; M.columns = [ALIAS.get(i, i) for i in M.columns]; M = M.loc[d.index, d.index]
    except FileNotFoundError: M = None
    return d, M
fig, ax = plt.subplots(2, 2, figsize=(11, 9.5)); rows = []
for c, (mol, proj) in enumerate((("DNA", a.dna), ("RNA", a.rna))):
    d, M = load(proj); rho, p = exact_spearman(d["within"].values, d["enriched_kmers"].values)
    rows.append(dict(molecule=mol, n_groups=len(d), spearman_rho=rho, exact_permutation_p=p)); t = d.assign(molecule=mol, ratio_between_within=d.mean_between / d.within); t.to_csv(f"{a.out}_{mol}_table.csv")
    if M is not None:
        im = ax[0, c].imshow(M.values, cmap="YlOrRd_r", vmin=M.values.min(), vmax=M.values.max()); ax[0, c].set_xticks(range(len(M))); ax[0, c].set_yticks(range(len(M)))
        ax[0, c].set_xticklabels(M.columns, rotation=60, ha="right"); ax[0, c].set_yticklabels(M.index)
        for i in range(len(M)):
            for j in range(len(M)): ax[0, c].text(j, i, f"{M.values[i, j]:.3f}", ha="center", va="center", fontsize=7.5, fontweight="bold" if i == j else "normal")
        ax[0, c].set_title(f"{mol}: mean Mash distance (k = {a.k}); diagonal = within group")
    else: ax[0, c].axis("off")
    ax[1, c].scatter(d["within"], d["enriched_kmers"], s=60, color="#2e75b6", zorder=3)
    for g, r in d.iterrows(): ax[1, c].annotate(g, (r["within"], r["enriched_kmers"]), textcoords="offset points", xytext=(5, 5), fontsize=9)
    ax[1, c].set_xlabel(f"Within-group Mash distance (k = {a.k})"); ax[1, c].set_ylabel("Number of enriched k-mers"); ax[1, c].grid(alpha=.3)
    ax[1, c].set_title(f"{mol}: Spearman \u03c1 = {rho:.2f}, exact permutation p = {p:.4f} (n = {len(d)})", fontsize=10)
for l, axx in zip("ABCD", ax.ravel()): axx.text(-0.12, 1.06, l, transform=axx.transAxes, fontsize=15, fontweight="bold")
plt.tight_layout(); plt.savefig(a.out + ".svg"); pd.DataFrame(rows).to_csv(a.out + "_spearman.csv", index=False); print(pd.DataFrame(rows).round(4).to_string(index=False))
