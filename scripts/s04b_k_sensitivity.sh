#!/usr/bin/env bash
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # folder of this script (works from any directory)
# usage: bash s04b_k_sensitivity.sh PROJECT_DIR k1 k2 ...   (needs the PATCHED portek on PATH)
P=$1; shift
for k in "$@"; do
  python "$SD/portek_patched_cli.py" find_k $P --min_k $k --max_k $k >/dev/null          # k-mer indices for this k only
  python "$SD/portek_patched_cli.py" find_enriched $P -k $k -m 2.0 > $P/output/log_k$k.txt
  python "$SD/portek_patched_cli.py" map $P -k $k -d 3 >> $P/output/log_k$k.txt
  cp $P/output/enriched_${k}mers_stats.csv $P/output/sens_enriched_k$k.csv
  cp $P/output/mapping_${k}mers_max_3_mismatches.tsv $P/output/sens_mapping_k$k.tsv
  echo "k=$k"; grep -E "enriched .*-mers" $P/output/log_k$k.txt
done
python - "$P" "$@" <<'PY'     # compare marker sets across k: counts per group and fraction mapped
import sys, pandas as pd
P, ks = sys.argv[1], sys.argv[2:]
rows = []
for k in ks:
    d = pd.read_csv(f"{P}/output/sens_enriched_k{k}.csv", index_col=0); m = pd.read_csv(f"{P}/output/sens_mapping_k{k}.tsv", sep="\t", index_col=0)
    e = d[d.group.str.endswith("_enriched")]; c = e.group.value_counts(); c["mapped_fraction"] = m.reference_sequence_position.notna().mean(); c.name = k; rows.append(c)
print(pd.DataFrame(rows).round(3))
PY
