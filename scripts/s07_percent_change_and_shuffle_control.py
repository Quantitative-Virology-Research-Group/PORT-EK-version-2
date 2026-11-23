#!/usr/bin/env python3
"""
(A) % change in enriched k-mers per group between two k values, from the s04b output files.
(B) chance-mapping control: map the real enriched k-mers and base-shuffled copies to the reference with the same index/d.
usage:
  python s07_percent_change_and_shuffle_control.py change PROJECT 15 17
  python s07_percent_change_and_shuffle_control.py shuffle PROJECT K D [n_kmers=2000] [n_shuffles=10] [seed=1]
(PORTEK_DIR must point to the PORT-EK-version-2 folder; the repository is not modified)
"""
import os, sys, random, yaml, numpy as np, pandas as pd
sys.path.insert(0, os.environ.get("PORTEK_DIR") or os.getcwd())
mode, proj = sys.argv[1], sys.argv[2]
if mode == "change":
    k1, k2 = sys.argv[3], sys.argv[4]
    cnt = {k: pd.read_csv(f"{proj}/output/sens_enriched_k{k}.csv", index_col=0).pipe(lambda d: d[d.group.str.endswith("_enriched")]).group.value_counts() for k in (k1, k2)}
    t = pd.DataFrame({f"k={k1}": cnt[k1], f"k={k2}": cnt[k2]}); t["pct_change"] = 100 * (t.iloc[:, 1] - t.iloc[:, 0]) / t.iloc[:, 0]
    print(t.round(1)); print(f"\nlargest absolute change k={k1}->{k2}: {t["pct_change"].abs().max():.1f}% ({t["pct_change"].abs().idxmax()})")
else:
    from Bio import SeqIO
    from portek.portek_map import MappingPipeline
    K, D = int(sys.argv[3]), int(sys.argv[4]); n = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
    ns = int(sys.argv[6]) if len(sys.argv) > 6 else 10; rnd = random.Random(int(sys.argv[7]) if len(sys.argv) > 7 else 1)
    cfg = yaml.safe_load(open(f"{proj}/config.yaml")); ref = str(next(SeqIO.parse(f"{proj}/input/{cfg['ref_seq']}", "fasta")).seq).upper()
    st = pd.read_csv(f"{proj}/output/enriched_{K}mers_stats.csv", index_col=0); st = st[st.group.str.endswith("_enriched")]
    kmers = list(st.sample(min(n, len(st)), random_state=1).index)
    mp = MappingPipeline.__new__(MappingPipeline); mp.k = K; mp.kmer_index = mp._generate_index(D, ref)
    def mapped(q): md = mp.map_kmer_to_index(q, D); return any(md[d] for d in range(D + 1))
    real = np.mean([mapped(q) for q in kmers]); sh = []
    for _ in range(ns):
        sh.append(np.mean([mapped("".join(rnd.sample(q, len(q)))) for q in kmers]))      # same base composition, random order
    print(f"k={K} d={D}: real mapped fraction = {real:.3f}; shuffled = {np.mean(sh):.4f} (SD {np.std(sh):.4f}, {ns} shuffles of {len(kmers)} k-mers)")
    print(f"estimated chance share of mapped k-mers = shuffled/real = {100*np.mean(sh)/real:.2f}%")
