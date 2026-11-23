#!/usr/bin/env python3
"""
R1-Major4 (isolate k-mer count) and R1-Major7 (random-walk / Markov analysis).
Input : <project>/output/<k>mer_counts_for_classifier.csv  and  enriched_<k>mers_stats.csv  (from `portek find_enriched`)
Usage : python s01_isolate_metric_and_graph.py PROJECT_DIR K [n_per_group=50] [n_perm=1000] [seed=42]
Replicates Analysis/K-mer.cout.space_Markov.chain.modeling.R (top-50 markers/group, exclusive first, kNN k=4,
edge weight = mean((x_i+x_j)/2), P = W/deg) and adds: identity check, stationary baseline, components,
label-permutation null, and the group-size-independent isolate metric.
"""
import sys, numpy as np, pandas as pd
from scipy.sparse.csgraph import connected_components
proj, K = sys.argv[1], int(sys.argv[2])
n_per = int(sys.argv[3]) if len(sys.argv) > 3 else 50
n_perm = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
rng = np.random.default_rng(int(sys.argv[5]) if len(sys.argv) > 5 else 42)
X = pd.read_csv(f"{proj}/output/{K}mer_counts_for_classifier.csv", index_col=0)
grp = X.pop("sample_group"); st = pd.read_csv(f"{proj}/output/enriched_{K}mers_stats.csv", index_col=0)
st = st[st.group.str.endswith("_enriched")]
groups = list(dict.fromkeys(grp))

# ---------- R1-Major4: isolate k-mer count ------------------------------------------------
N = grp.map(grp.value_counts())
orig = X.sum(axis=1) / N                      # what the R code does: rowSums(all enriched k-mers)/N_group
raw  = X.sum(axis=1)                          # group-size independent
own = pd.Series({s: (X.loc[s, st.index[st.group == f"{grp[s]}_enriched"].intersection(X.columns)] > 0).mean()
                 for s in X.index})           # fraction of the isolate's OWN-group markers present
fracall = (X > 0).mean(axis=1)               # fraction of all enriched k-mers present
out = pd.DataFrame({"group": grp, "N_group": N, "orig_sum_over_N": orig, "raw_sum": raw,
                    "frac_own_group_markers": own, "frac_all_markers": fracall})
print("== Isolate metric by group (median) ==")
print(out.groupby("group").median(numeric_only=True).round(3))
out.to_csv(f"{proj}/output/isolate_metric_variants.csv")

# ---------- R1-Major7: graph + Markov chain ------------------------------------------------
def top_markers(g, n=50):
    c = st[st.group == f"{g}_enriched"].copy()
    c["excl"] = (c.exclusivity == "exclusive").astype(int)
    return list(c.sort_values(["excl", f"{g}_avg", "RMSE"], ascending=[False, False, True]).head(n).index)
markers = list(dict.fromkeys(m for g in groups for m in top_markers(g)))
markers = [m for m in markers if m in X.columns]
idx = np.concatenate([rng.choice(np.where(grp == g)[0], min(n_per, (grp == g).sum()), replace=False) for g in groups])
F = X.iloc[idx][markers].to_numpy(float); lab = grp.iloc[idx].to_numpy()
n = len(F); d = np.sqrt(((F[:, None, :] - F[None, :, :]) ** 2).sum(-1)); np.fill_diagonal(d, np.inf)
W = np.zeros((n, n))
for i in range(n):
    for j in np.argsort(d[i])[:4]:                    # k = 4 nearest neighbours (undirected union)
        w = ((F[i] + F[j]) / 2).mean(); W[i, j] = W[j, i] = w
deg = W.sum(1); deg[deg == 0] = 1e-9; P = W / deg[:, None]
def stats(labels):
    same = labels[:, None] == labels[None, :]; np.fill_diagonal(same, False)
    pw = (P * same).sum(1); return pw
pw = stats(lab); print("\n== Random-walk summary ==")
print(f"markers={len(markers)}  nodes={n}  mean p_within={pw.mean():.3f}  p_across={1-pw.mean():.3f}")
# (i) algebraic identity: with equal n/group, mean(diag of lumped matrix) == mean(p_within)
Pl = np.array([[P[np.ix_(lab == s, lab == t)].sum(1).mean() for t in groups] for s in groups])
print("group order:", groups)
print("lumped diag (same order):", dict(zip(groups, np.round(np.diag(Pl), 3))), " mean diag =", round(np.diag(Pl).mean(), 3),
      "(== mean p_within when groups have equal n; any difference in the paper = different run/aggregation)")
# (ii) connectivity
nc, comp = connected_components(W > 0); sz = np.bincount(comp); print("connected components:", nc, "sizes:", sz)
iso = (W.sum(1) == 0); print(f"isolated nodes (no edge weight, e.g. isolates carrying none of the markers): {int(iso.sum())}; mean p_within excluding them = {pw[~iso].mean():.3f}")
for c_ in np.where(sz < 0.5 * n)[0]: print(f"  component {c_}: size {sz[c_]}, groups: {dict(zip(*np.unique(lab[comp == c_], return_counts=True)))}")
# (iii) correct baseline: stationary mass of each group, pi_i = deg_i / sum(deg)  (not 1/7)
pi = deg / deg.sum(); base = np.mean([pi[lab == g].sum() for g in groups]); print("stationary (degree-proportional) mass per group:", {g: round(float(pi[lab == g].sum()), 3) for g in groups}, "(meaningful only within connected components)")
# (iv) label-permutation null on the SAME graph (permute group labels over nodes)
null = np.array([stats(rng.permutation(lab)).mean() for _ in range(n_perm)])
print(f"permutation null mean p_within = {null.mean():.3f} (sd {null.std():.3f}); observed = {pw.mean():.3f}; "
      f"empirical p = {(1 + (null >= pw.mean()).sum()) / (1 + n_perm):.4f}")
# (v) multi-step: exact occupancy of home group via matrix powers vs 1/7 and vs stationary baseline
Mhome = (lab[:, None] == lab[None, :]).astype(float); Pt = np.eye(n); occ = []
for t in range(50): Pt = Pt @ P; occ.append((Pt * Mhome).sum(1))
obs_occ = np.mean(occ); print(f"50-step time-averaged home occupancy = {obs_occ:.3f}  (uniform 1/7 = {1/len(groups):.3f}; the mean stationary mass is always 1/7 and is not a useful baseline)")
Ps = []; Pt = np.eye(n)
for t_ in range(50): Pt = Pt @ P; Ps.append(Pt.copy())
nulls = []
for _ in range(min(n_perm, 200)):
    pl = rng.permutation(lab); Mh = (pl[:, None] == pl[None, :]).astype(float); nulls.append(np.mean([(Q * Mh).sum(1).mean() for Q in Ps]))
nulls = np.array(nulls); print(f"permutation null of 50-step occupancy ({len(nulls)} label shuffles on the same graph): mean = {nulls.mean():.3f} (sd {nulls.std():.3f}); observed = {obs_occ:.3f}; empirical p = {(1 + (nulls >= obs_occ).sum()) / (1 + len(nulls)):.4f}")
# (vi) does lumping reproduce t-step dynamics?  compare Pl^t with exact group-level t-step averages
for t in (2, 5, 10):
    Ex = np.linalg.matrix_power(P, t); true_t = np.array([[Ex[np.ix_(lab == s, lab == u)].sum(1).mean() for u in groups] for s in groups])
    print(f"t={t}: max|Pl^t - exact| = {np.abs(np.linalg.matrix_power(Pl, t) - true_t).max():.3f}")
