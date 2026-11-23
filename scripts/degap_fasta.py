#!/usr/bin/env python3
"""Remove alignment gaps ('-') from LANL FASTA files so that PORT-EK-v2 computes L (find_k) from real nucleotides.
usage: python degap_fasta.py IN_DIR OUT_DIR  (all *.fasta)"""
import sys, glob, os
from Bio import SeqIO
os.makedirs(sys.argv[2], exist_ok=True)
for f in glob.glob(f"{sys.argv[1]}/*.fasta"):
    recs = list(SeqIO.parse(f, "fasta"))
    for r in recs: r.seq = r.seq.upper().replace("-", ""); r.description = r.id
    SeqIO.write(recs, f"{sys.argv[2]}/{os.path.basename(f)}", "fasta"); print(os.path.basename(f), len(recs))
