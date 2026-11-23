#!/usr/bin/env python3
"""Create a balanced random-subsample PORT-EK project (for sensitivity analyses; R1 Major 3).
usage: python scripts/setup_sub.py DNA 80 --seed 1 [--data_dir data/DNA_degapped] [--out_dir PROJ] [--ref HXB2.fasta]
  mol      DNA or RNA            n   sequences per group (all if a group has fewer)
  --seed   an integer (1, 2, 3 ...); loop over seeds for replicates:  for s in $(seq 1 100); do python scripts/setup_sub.py DNA 80 --seed $s; done
Input files must be named hiv-db.<DNA|RNA>.<A.comp|B|C|D|Rare|CRF01.AE|CRF02.AG>.fasta (as in the repository zips). Gaps are removed."""
import argparse, os, random, yaml
from Bio import SeqIO
ap = argparse.ArgumentParser(); ap.add_argument("mol", choices=["DNA", "RNA"]); ap.add_argument("n", type=int)
ap.add_argument("--seed", type=int, default=1); ap.add_argument("--data_dir"); ap.add_argument("--out_dir"); ap.add_argument("--dedup", action="store_true", help="remove identical (ungapped) sequences before sampling"); ap.add_argument("--ref", help="HXB2 FASTA (K03455); default: derived from the B file")
a = ap.parse_args(); mol = a.mol; random.seed(a.seed)
data = a.data_dir or f"data/{mol}_degapped"; proj = a.out_dir or f"proj_{mol}_{a.n}_seed{a.seed}"
if not os.path.isdir(data): raise SystemExit(f"data directory not found: {data} (use --data_dir)")
for s in ("input", "output", "temp"): os.makedirs(f"{proj}/{s}", exist_ok=True)
groups = ["A", "B", "C", "D", "Rare", "CRF01.AE", "CRF02.AG"]; fm = {"A": "A.comp", "B": "B", "C": "C", "D": "D", "Rare": "Rare", "CRF01.AE": "CRF01.AE", "CRF02.AG": "CRF02.AG"}
files = []
for g in groups:
    p = f"{data}/hiv-db.{mol}.{fm[g]}.fasta"
    if not os.path.exists(p): raise SystemExit(f"missing file: {p}")
    recs = list(SeqIO.parse(p, "fasta"))
    if a.dedup:
        seen = set(); keep = []
        for r in recs:
            h = str(r.seq).upper().replace("-", "")
            if h not in seen: seen.add(h); keep.append(r)
        print(g, "deduplicated:", len(recs), "->", len(keep)); recs = keep
    sel = random.sample(recs, min(a.n, len(recs)))
    for r in sel: r.seq = r.seq.upper().replace("-", ""); r.description = r.id
    SeqIO.write(sel, f"{proj}/input/{g}.fasta", "fasta"); files.append(f"{g}.fasta"); print(g, len(recs), "->", len(sel))
if a.ref: open(f"{proj}/input/HXB2.fasta", "w").write(open(a.ref).read())
else:
    for r in SeqIO.parse(f"{data}/hiv-db.{mol}.B.fasta", "fasta"):
        if "HXB2" in r.id or "K03455" in r.id:
            open(f"{proj}/input/HXB2.fasta", "w").write(f">HXB2\n{str(r.seq).upper().replace('-', '')}\n"); print("WARNING: reference taken from", r.id, "- pass --ref with K03455 for the paper analysis"); break
yaml.safe_dump({"sample_groups": groups, "input_files": files, "header_format": [], "mode": "ava", "goi": None, "ref_seq": "HXB2.fasta"}, open(f"{proj}/config.yaml", "w"))
print("project written to", proj)
