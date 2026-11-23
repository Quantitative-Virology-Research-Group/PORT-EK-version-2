#!/usr/bin/env python3
"""
R1-Major8: independent within-/between-group diversity (alignment-free Mash-type distance) to relate to marker counts.
Distance:  D = -ln(2J/(1+J))/k  with J = Jaccard of the two sequences' k-mer sets (Ondov et al. 2016); ~ per-site divergence.
Usage: python s03_diversity_kmer_mash.py PROJECT_DIR [k=13] [pairs=1500] [seed=1]
(Use k=11-13 for between-subtype divergence; k=17 saturates. For the gold-standard, also run MAFFT --auto per group
 and compute mean pairwise p-distance, see the Word document, step 8.)
"""
import sys, itertools, numpy as np, yaml
from Bio import SeqIO
proj = sys.argv[1]; k = int(sys.argv[2]) if len(sys.argv) > 2 else 13
npairs = int(sys.argv[3]) if len(sys.argv) > 3 else 1500; rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
cfg = yaml.safe_load(open(f"{proj}/config.yaml")); G = cfg["sample_groups"]
code = np.full(256, 255, np.uint8)
for i, c in enumerate("ACGT"): code[ord(c)] = i
def kmerset(seq):
    a = code[np.frombuffer(str(seq).upper().replace("-", "X").encode(), np.uint8)]          # gaps/ambiguity break k-mers (as PORT-EK)
    n = len(a) - k + 1
    if n <= 0: return np.array([], np.uint64)
    ok = np.ones(n, bool); h = np.zeros(n, np.uint64)
    for j in range(k):
        s = a[j:j + n]; ok &= s < 4; h = (h << np.uint64(2)) | (s & 3).astype(np.uint64)
    return np.unique(h[ok])
S = {g: [kmerset(r.seq) for r in SeqIO.parse(f"{proj}/input/{f}", "fasta")] for g, f in zip(G, cfg["input_files"])}
def mash(a, b):
    inter = np.intersect1d(a, b, assume_unique=True).size; J = inter / (a.size + b.size - inter)
    return np.nan if J == 0 else -np.log(2 * J / (1 + J)) / k
def mean_d(A, B, same):
    if same: pairs = list(itertools.combinations(range(len(A)), 2))
    else: pairs = [(i, j) for i in range(len(A)) for j in range(len(B))]
    if len(pairs) > npairs: pairs = [pairs[i] for i in rng.choice(len(pairs), npairs, replace=False)]
    return np.nanmean([mash(A[i], B[j]) for i, j in pairs])
import pandas as pd
M = pd.DataFrame(index=G, columns=G, dtype=float)
for a in G:
    for b in G:
        if G.index(b) >= G.index(a): M.loc[a, b] = M.loc[b, a] = mean_d(S[a], S[b], a == b)
print(f"Mean Mash distance (k={k}); diagonal = within-group diversity, off-diagonal = between-group"); print(M.round(3))
within = pd.Series(np.diag(M), index=G); between = (M.sum(axis=1) - within) / (len(G) - 1)
M.to_csv(f"{proj}/output/diversity_mash_matrix_k{k}.csv"); res = pd.DataFrame({"within": within, "mean_between": between}); res.to_csv(f"{proj}/output/diversity_mash_k{k}.csv"); print(res.round(3))
