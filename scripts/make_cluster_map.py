#!/usr/bin/env python3
"""Convert a cd-hit-est .clstr file to a 2-column map (sequence_id<TAB>cluster_id).
usage: python make_cluster_map.py clusters.clstr project_rebuttal/cluster_map.tsv
Run cd-hit-est with -d 0, otherwise IDs longer than ~19 characters are truncated in the .clstr:
  cd-hit-est -i all_degapped.fasta -o clusters -c 0.98 -n 10 -d 0 -T 8 -M 0"""
import re, sys
cl = None; n = 0
with open(sys.argv[1]) as f, open(sys.argv[2], "w") as o:
    for l in f:
        if l.startswith(">Cluster"): cl = l.split()[1]
        else:
            m = re.search(r">(.+?)\.\.\.( |$)", l)
            if m: o.write(f"{m.group(1)}\t{cl}\n"); n += 1
print(n, "sequences written to", sys.argv[2])
