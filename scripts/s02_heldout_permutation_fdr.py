#!/usr/bin/env python3
"""
R1-Major1/2/3 + R2-4: held-out validation, label-permutation null of the WHOLE selection procedure,
BY-FDR vs raw p-values, effect size and bootstrap stability.
Two stages (so each fits in a normal session):
  python s02_heldout_permutation_fdr.py prepare PROJECT -k 15 --cluster_map DNA_degapped_clusters.tsv [--n_disc 0] [--max_mem 2.0] [--n_jobs 8]
  python s02_heldout_permutation_fdr.py analyze PROJECT_discall -k 15 --n_perm 100 --n_boot 100
Uses the ORIGINAL portek + portek_runtime_patch.py (in-memory fix; no repository file is edited). Split unit here = sequence; for the paper split by
patient/study (LANL metadata) or by sequence cluster (see Word document) and pass the resulting FASTA lists instead.
"""
import sys, os, random, shutil, pickle, yaml, numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.environ.get("PORTEK_DIR") or os.getcwd())   # ORIGINAL repository, unmodified
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from Bio import SeqIO
import portek, portek_runtime_patch   # in-memory fix of the group assignment (repository files untouched)

def load_cluster_map(path, proj):
    """Accepts (a) 2-column TSV  id<TAB>cluster  or (b) a cd-hit-est .clstr file. Path may be relative to the project folder."""
    import re
    p = path if os.path.exists(path) else os.path.join(proj, path)
    if not os.path.exists(p): sys.exit(f"cluster map not found: tried '{path}' and '{p}'")
    first = open(p).readline()
    if first.startswith(">") and not first.startswith(">Cluster") and os.path.exists(p + ".clstr"):   # cd-hit-est -o X writes X (FASTA) and X.clstr
        print(f"NOTE: {p} is cd-hit's FASTA output; using {p}.clstr instead"); p += ".clstr"; first = open(p).readline()
    if first.startswith(">Cluster"):                                   # cd-hit .clstr
        parent = {}
        def find(x):
            while parent.setdefault(x, x) != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        cm, cl, merged = {}, None, 0
        for l in open(p):
            if l.startswith(">Cluster"): cl = l.split()[1]; find(cl)
            else:
                m = re.search(r">(.+?)\.\.\.( |$)", l)
                if m:
                    k = m.group(1)
                    if k in cm and find(cm[k]) != find(cl): parent[find(cl)] = find(cm[k]); merged += 1   # identical (e.g. truncated) IDs in different clusters -> merge (conservative)
                    cm.setdefault(k, cl)
        if merged: print(f"NOTE: {merged} cluster merges because several sequences share the same ID string in the .clstr (truncated IDs?)")
        return {k: find(c) for k, c in cm.items()}
    if first.startswith(">"): sys.exit(f"{p} is a FASTA file (cd-hit-est's representative sequences), not a cluster map. Use the '.clstr' file written next to it (e.g. clusters.clstr), and run cd-hit-est with -d 0 so full IDs are kept.")
    cm = {}
    for l in open(p):
        a = l.rstrip("\n").split("\t")
        if len(a) >= 2: cm[a[0]] = a[1]
    if not cm: sys.exit(f"no 'id<TAB>cluster' lines found in {p}")
    return cm

def reconcile_ids(cm, ids):
    """If cd-hit truncated the IDs (run without -d 0), match records by the truncated prefix; ambiguous prefixes were already merged."""
    miss = [i for i in ids if i not in cm]
    if not miss: return cm
    L = max(len(k) for k in cm)
    if L > 20: sys.exit(f"{len(miss)} sequence IDs are not in the cluster map, e.g. {miss[:3]}; map IDs are up to {L} characters, so this is not truncation. IDs must match the FASTA headers exactly.")
    out, bad = {}, []
    for i in ids:
        k = i if i in cm else i[:L]
        if k in cm: out[i] = cm[k]
        else: bad.append(i)
    if bad: sys.exit(f"{len(bad)} IDs cannot be matched even by {L}-character prefix, e.g. {bad[:3]}")
    print(f"WARNING: cd-hit truncated IDs to {L} characters; matched {len(miss)} sequences by prefix (clusters sharing a prefix were merged, which only makes the split stricter). Rerun cd-hit-est with -d 0 for exact matching.")
    return out

VERSION = "s02 v2026-11-03 (cluster .clstr + truncated-ID fallback)"

def prepare(src, K, n_disc=0, seed=1, cluster_map=None, max_mem=2.0, min_freq=0.1, n_jobs=1):
    print(VERSION, flush=True)
    if not os.path.exists(f"{src}/config.yaml"): sys.exit(f"{src}/config.yaml not found -- PROJECT must be the folder with config.yaml and input/")
    cfg = yaml.safe_load(open(f"{src}/config.yaml")); random.seed(seed)
    disc, test = f"{src}_disc{n_disc or 'all'}", f"{src}_test"
    for d in (disc, test):
        shutil.rmtree(d, ignore_errors=True)
        for s in ("input", "output", "temp"): os.makedirs(f"{d}/{s}")
    shutil.copy(f"{src}/input/{cfg['ref_seq']}", f"{disc}/input/"); yaml.safe_dump(cfg, open(f"{disc}/config.yaml", "w"))
    cm_all = load_cluster_map(cluster_map, src) if cluster_map else None
    if cm_all: cm_all = reconcile_ids(cm_all, [r.id for f_ in cfg["input_files"] for r in SeqIO.parse(f"{src}/input/{f_}", "fasta")])
    for g, f in zip(cfg["sample_groups"], cfg["input_files"]):
        recs = list(SeqIO.parse(f"{src}/input/{f}", "fasta")); random.shuffle(recs)
        if cluster_map:                         # whole clusters go to the held-out set (25% of clusters)
            cm = cm_all
            miss = [r.id for r in recs if r.id not in cm]
            if miss:
                L = max(len(k) for k in cm); hint = f" The longest ID in the map has {L} characters: cd-hit truncated the IDs -- rerun cd-hit-est with -d 0." if L <= 20 else ""
                sys.exit(f"{len(miss)} sequence IDs of {f} are not in the cluster map, e.g. {miss[:3]}.{hint} IDs must match the FASTA headers exactly.")
            cl = sorted({cm[r.id] for r in recs}); random.Random(seed).shuffle(cl)
            tcl = set(cl[:max(1, int(0.25 * len(cl)))]); tr = [r for r in recs if cm[r.id] not in tcl]; tr = tr[:n_disc] if n_disc > 0 else tr; te = [r for r in recs if cm[r.id] in tcl]
        else:
            if n_disc <= 0: sys.exit('random split needs --n_disc > 0 (or give --cluster_map)')
            tr, te = recs[:n_disc], recs[n_disc:]
        SeqIO.write(tr, f"{disc}/input/{f}", "fasta"); SeqIO.write(te, f"{test}/input/{f}", "fasta")
        print(g, "discovery", len(tr), "held-out", len(te))
    portek.KmerFinder(disc, K, K).find_all_kmers(n_jobs=n_jobs)
    e = portek.EnrichedKmersPipeline(disc, K); e.get_basic_kmer_stats(max_mem=max_mem, min_freq=min_freq)
    e.calc_kmer_stats("common", n_jobs=n_jobs, false_discovery_control=False)
    m = e.matrices["common"]; samples = [c for c in m.columns if c in set(e.sample_list)]
    C = m[samples].to_numpy(np.uint8); lab = np.array([s.split("_")[0] for s in samples])
    other = m.drop(columns=samples)
    pickle.dump(dict(C=C, lab=lab, idx=m.index.to_numpy(), groups=cfg["sample_groups"], stats=other), open(f"{disc}/output/common_matrix.pkl", "wb"))
    e.get_enriched_kmers(); print("pipeline enriched:\n", e.matrices["enriched"].group.value_counts())
    print(f"\nDONE. Discovery project: {disc}   Held-out sequences: {test}\nNext: python scripts/s02_heldout_permutation_fdr.py analyze {disc} -k {K}", flush=True)

def enrich(C, lab, groups, fdr=False, rmse_min=0.1, pthr=0.01):
    """Vectorised re-implementation of PORT-EK-v2 'ava' selection (with exact group matching). Returns assigned group index or -1."""
    gi = [np.where(lab == g)[0] for g in groups]; s = len(groups)
    avg = np.stack([C[:, i].mean(1) for i in gi], 1)
    pairs = [(a, b) for b in range(1, s) for a in range(b)]
    P = np.empty((C.shape[0], len(pairs))); E = np.empty_like(P)
    for j, (a, b) in enumerate(pairs):
        P[:, j] = stats.mannwhitneyu(C[:, gi[a]], C[:, gi[b]], axis=1).pvalue; E[:, j] = avg[:, a] - avg[:, b]
        if fdr: P[:, j] = stats.false_discovery_control(P[:, j], method="by")
    rmse = np.sqrt((E ** 2).mean(1)); cand = avg.argmax(1)
    inv = np.array([[g in pr for pr in pairs] for g in range(s)])            # exact membership (the fix)
    ok = np.array([(P[r, inv[cand[r]]] < pthr).all() for r in range(C.shape[0])]) & (rmse > rmse_min)
    return np.where(ok, cand, -1), avg, P, E, rmse, pairs

def analyze(disc, K, n_perm, n_boot, sub=20000):
    print(VERSION, flush=True)
    if not os.path.exists(f"{disc}/output/common_matrix.pkl"):
        import glob; found = sorted(glob.glob(os.path.dirname(os.path.abspath(disc)) + "/*_disc*"))
        sys.exit(f"{disc}/output/common_matrix.pkl not found: `prepare` has not finished successfully (it prints DONE at the end). Existing discovery folders: {found or 'none'}. Rerun prepare and read its last lines.")
    d = pickle.load(open(f"{disc}/output/common_matrix.pkl", "rb")); C, lab, groups = d["C"], d["lab"], d["groups"]
    rng = np.random.default_rng(0); print(f"tested k-mers (after entropy/min_freq filter) = {C.shape[0]}; samples = {C.shape[1]}")
    # 1) raw vs BY-FDR
    a_raw, avg, P, E, rmse, pairs = enrich(C, lab, groups, False); a_fdr = enrich(C, lab, groups, True)[0]
    t = pd.DataFrame({"raw_p<0.01": pd.Series(a_raw[a_raw >= 0]).map(dict(enumerate(groups))).value_counts(),
                      "BY-FDR<0.01": pd.Series(a_fdr[a_fdr >= 0]).map(dict(enumerate(groups))).value_counts()}).fillna(0).astype(int)
    print("\nEnriched k-mers per group: raw vs Benjamini-Yekutieli\n", t, "\nretained after FDR: %.1f%%" % (100 * (a_fdr >= 0).sum() / max(1, (a_raw >= 0).sum())))
    t.to_csv(f"{disc}/output/enriched_raw_vs_fdr.csv")
    # 2) label-permutation null of the full procedure (on a random k-mer subset; scale by C.shape[0]/sub)
    sel = rng.choice(C.shape[0], min(sub, C.shape[0]), replace=False); scale = C.shape[0] / len(sel)
    obs_sub = (enrich(C[sel], lab, groups)[0] >= 0).sum() * scale
    null = [(enrich(C[sel], rng.permutation(lab), groups)[0] >= 0).sum() * scale for _ in range(n_perm)]
    print(f"\nPermutation null (whole selection, {n_perm} shuffles, k-mer-subsample scaled): mean={np.mean(null):.0f} max={np.max(null):.0f}; observed={obs_sub:.0f}")
    # 3) bootstrap stability of observed markers (subset of 2000 markers)
    mk = rng.choice(np.where(a_raw >= 0)[0], min(2000, (a_raw >= 0).sum()), replace=False); keep = np.zeros(len(mk)); 
    for b in range(n_boot):
        bs = np.concatenate([rng.choice(np.where(lab == g)[0], (lab == g).sum()) for g in groups])
        keep += (enrich(C[mk][:, bs], lab[bs], groups)[0] == a_raw[mk])
    stab = keep / n_boot; print(f"\nBootstrap stability (n={n_boot}) of 2000 random markers: median={np.median(stab):.2f}; fraction >=0.8: {(stab >= 0.8).mean():.2f}")
    # 4) effect size: Cliff's delta (candidate vs mean of other groups) for markers
    # 5) held-out classification with frozen markers
    markers = {g: d["idx"][a_raw == i] for i, g in enumerate(groups)}
    Xd = {int(m): (g) for g, ms in markers.items() for m in ms}
    test = f"{disc.split('_disc')[0]}_test"; cfg = yaml.safe_load(open(f"{disc}/config.yaml")); K_ = K
    code = {c: i for i, c in enumerate("ACGT")}; y, pred = [], []
    mset = {int(m): g for g, ms in markers.items() for m in ms}; msize = {g: len(ms) for g, ms in markers.items()}
    for g, f in zip(groups, cfg["input_files"]):
        for r in SeqIO.parse(f"{test}/input/{f}", "fasta"):
            s = str(r.seq).upper(); score = {gg: 0 for gg in groups}; seen = set()
            for i in range(len(s) - K_ + 1):
                w = s[i:i + K_]
                if all(c in code for c in w):
                    v = 0
                    for c in w: v = (v << 2) | code[c]
                    if v in mset and v not in seen: seen.add(v); score[mset[v]] += 1
            norm = {gg: score[gg] / max(1, msize[gg]) for gg in groups}; pred.append(max(norm, key=norm.get)); y.append(g)
    from sklearn.metrics import confusion_matrix, classification_report, f1_score
    print("\nHeld-out assignment (rule: group with highest fraction of its own markers present)")
    print(pd.DataFrame(confusion_matrix(y, pred, labels=groups), index=groups, columns=groups))
    print(classification_report(y, pred, labels=groups, digits=3, zero_division=0)); print("macro-F1 = %.3f" % f1_score(y, pred, average="macro", labels=groups))

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Held-out validation, permutation null, BY-FDR, bootstrap (original portek + runtime patch)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("prepare", help="split sequences and run discovery"); a.add_argument("project"); a.add_argument("-k", type=int, required=True)
    a.add_argument("--n_disc", type=int, default=0, help="max discovery sequences per group (0 = all non-test; required >0 for a random split)")
    a.add_argument("--cluster_map", help="TSV: sequence_id<TAB>cluster; whole clusters go to the held-out set (25%%)"); a.add_argument("--seed", type=int, default=1)
    a.add_argument("--max_mem", type=float, default=2.0); a.add_argument("--min_freq", type=float, default=0.1); a.add_argument("--n_jobs", type=int, default=1)
    b = sub.add_parser("analyze", help="FDR, permutation null, bootstrap, held-out scoring"); b.add_argument("disc_project"); b.add_argument("-k", type=int, required=True)
    b.add_argument("--n_perm", type=int, default=100); b.add_argument("--n_boot", type=int, default=100); b.add_argument("--kmer_subsample", type=int, default=20000)
    r = ap.parse_args()
    if r.cmd == "prepare": prepare(r.project, r.k, r.n_disc, r.seed, r.cluster_map, r.max_mem, r.min_freq, r.n_jobs)
    else: analyze(r.disc_project, r.k, r.n_perm, r.n_boot, r.kmer_subsample)
