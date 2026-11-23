#!/usr/bin/env python3
"""
R1-Major3: balanced downsampling replicates (deduplicated) -> corrected PORT-EK-v2 enrichment on each replicate -> summary.
Runs the ORIGINAL pipeline in-process (no plots) with portek_runtime_patch.py; repository files are not edited.

usage (from the PORT-EK-version-2 folder):
  python scripts/s08_replicates.py DNA -k 15 --n_per_group 78 --replicates 100 --workers 4 \
      --data_dir data/DNA_degapped --full_stats project_rebuttal/output/enriched_15mers_stats.csv --max_mem 2.0 --out reps_DNA_k15
--n_per_group 0  -> use the size of the smallest group after deduplication (balanced exactly)
--no_dedup       -> skip removal of identical (ungapped) sequences
Memory: about 1-2 GB per worker for 7 x 80 genomes.  Restartable: finished replicates (out/rep_XXX.pkl) are skipped.
"""
import argparse, os, sys, shutil, random, pickle, hashlib, yaml, glob, time
from concurrent.futures import ProcessPoolExecutor, as_completed
GROUPS = ["A", "B", "C", "D", "Rare", "CRF01.AE", "CRF02.AG"]; FM = {"A": "A.comp", "B": "B", "C": "C", "D": "D", "Rare": "Rare", "CRF01.AE": "CRF01.AE", "CRF02.AG": "CRF02.AG"}

def load_groups(mol, data_dir, dedup):
    from Bio import SeqIO
    out, sizes = {}, {}
    for g in GROUPS:
        recs, seen = [], set()
        for r in SeqIO.parse(f"{data_dir}/hiv-db.{mol}.{FM[g]}.fasta", "fasta"):
            s = str(r.seq).upper().replace("-", "")
            if dedup:
                h = hashlib.md5(s.encode()).hexdigest()
                if h in seen: continue
                seen.add(h)
            recs.append((r.id, s))
        out[g] = recs; sizes[g] = len(recs)
    return out, sizes

def one_replicate(args):
    seed, mol, k, n, data_dir, dedup, out, max_mem, min_freq, keep = args
    res = f"{out}/rep_{seed:03d}.pkl"
    if os.path.exists(res): return seed, "skipped"
    sys.path.insert(0, os.environ.get("PORTEK_DIR") or os.getcwd()); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import portek, portek_runtime_patch, pandas as pd
    groups, sizes = load_groups(mol, data_dir, dedup); n_eff = n if n > 0 else min(sizes.values())
    rnd = random.Random(seed); proj = f"{out}/work_{seed:03d}"
    for s in ("input", "output", "temp"): os.makedirs(f"{proj}/{s}", exist_ok=True)
    files = []
    for g in GROUPS:
        sel = rnd.sample(groups[g], min(n_eff, len(groups[g])))
        with open(f"{proj}/input/{g}.fasta", "w") as f:
            for i, s in sel: f.write(f">{i}\n{s}\n")
        files.append(f"{g}.fasta")
    open(f"{proj}/input/HXB2.fasta", "w").write(">ref\nACGT\n")
    yaml.safe_dump({"sample_groups": GROUPS, "input_files": files, "header_format": [], "mode": "ava", "goi": None, "ref_seq": "HXB2.fasta"}, open(f"{proj}/config.yaml", "w"))
    portek.KmerFinder(proj, k, k).find_all_kmers(n_jobs=1)
    e = portek.EnrichedKmersPipeline(proj, k); e.get_basic_kmer_stats(max_mem=max_mem, min_freq=min_freq)
    e.calc_kmer_stats("common", n_jobs=1, false_discovery_control=False); tested = len(e.matrices["common"])
    e.get_enriched_kmers(); m = e.matrices["enriched"]; m = m[m.group.str.endswith("_enriched")]
    pickle.dump({"seed": seed, "n_per_group": n_eff, "tested": tested, "group": {portek.decode_kmer(i, k): g for i, g in m.group.str.replace("_enriched", "").items()}}, open(res, "wb"))
    if not keep: shutil.rmtree(proj, ignore_errors=True)
    return seed, f"done tested={tested} enriched={len(m)}"

def decode_kmer(i, k):
    s = bin(int(i))[2:].rjust(2 * k, "0"); return "".join({"00": "A", "01": "C", "10": "G", "11": "T"}[s[j:j + 2]] for j in range(0, 2 * k, 2))

def summarise(out, k, full_stats):
    import numpy as np, pandas as pd
    reps = [pickle.load(open(f, "rb")) for f in sorted(glob.glob(f"{out}/rep_*.pkl"))]
    if not reps: return
    for r in reps:     # replicates saved by earlier versions hold integer-coded k-mers; the full-data stats file holds sequences
        r["group"] = {(decode_kmer(km, k) if not isinstance(km, str) else km): g for km, g in r["group"].items()}
    cnt = pd.DataFrame([pd.Series(r["group"]).value_counts() for r in reps]).reindex(columns=GROUPS).fillna(0).astype(int)
    cnt["tested"] = [r["tested"] for r in reps]; cnt.index = [r["seed"] for r in reps]; cnt.to_csv(f"{out}/replicate_counts.csv")
    d = cnt.describe(percentiles=[.025, .975]).loc[["mean", "std", "2.5%", "97.5%", "min", "max"]].T
    print(f"\n{len(reps)} replicates, n per group = {reps[0]['n_per_group']} (after deduplication)\nEnriched k-mers per group (mean, SD, 2.5-97.5 percentile):\n", d.round(1))
    d.to_csv(f"{out}/replicate_summary.csv")
    allk = {}
    for r in reps:
        for km, g in r["group"].items(): allk.setdefault((km, g), 0); allk[(km, g)] += 1
    rec = pd.Series(allk); rec.index = pd.MultiIndex.from_tuples(rec.index, names=["kmer", "group"]); rec = rec / len(reps)
    print("\nRecurrence of k-mers across replicates (per group: union size; share enriched in >=50% / >=90% of replicates):")
    print(rec.groupby("group").agg(union="size", ge50=lambda x: (x >= .5).mean(), ge90=lambda x: (x >= .9).mean()).round(3))
    rec.to_csv(f"{out}/kmer_recurrence.csv")
    if full_stats and os.path.exists(full_stats):
        f = pd.read_csv(full_stats, index_col=0); f = f[f.group.str.endswith("_enriched")]; fg = f.group.str.replace("_enriched", "").replace({"R": "Rare"})   # full run may name the Rare pool "R"
        rows = []
        for i, r in enumerate(reps):
            s = pd.Series(r["group"]); ok = s.index.isin(fg.index) & (fg.reindex(s.index).values == s.values)
            rows.append(pd.Series(ok, index=s.index).groupby(s.values).mean())
        c = pd.DataFrame(rows).reindex(columns=GROUPS)
        print("\nShare of replicate markers also enriched for the SAME group in the full-data run (mean over replicates):\n", c.mean().round(3)); c.to_csv(f"{out}/replicate_vs_full.csv")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("mol", choices=["DNA", "RNA"]); ap.add_argument("-k", type=int, required=True)
    ap.add_argument("--n_per_group", type=int, default=0); ap.add_argument("--replicates", type=int, default=100); ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--data_dir"); ap.add_argument("--full_stats"); ap.add_argument("--max_mem", type=float, default=2.0); ap.add_argument("--min_freq", type=float, default=0.1)
    ap.add_argument("--no_dedup", action="store_true"); ap.add_argument("--keep", action="store_true"); ap.add_argument("--out"); ap.add_argument("--summary_only", action="store_true")
    a = ap.parse_args(); data = a.data_dir or f"data/{a.mol}_degapped"; out = a.out or f"reps_{a.mol}_k{a.k}"; os.makedirs(out, exist_ok=True)
    if not a.summary_only:
        g, sizes = load_groups(a.mol, data, not a.no_dedup); print("group sizes after", "no" if a.no_dedup else "", "deduplication:", sizes, flush=True)
        jobs = [(s, a.mol, a.k, a.n_per_group, data, not a.no_dedup, out, a.max_mem, a.min_freq, a.keep) for s in range(1, a.replicates + 1)]
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for fu in as_completed([ex.submit(one_replicate, j) for j in jobs]): print(f"[{time.time()-t0:6.0f}s] replicate", *fu.result(), flush=True)
    summarise(out, a.k, a.full_stats)
