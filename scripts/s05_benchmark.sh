#!/usr/bin/env bash
# R1-Major5 / R2-1..3: like-for-like benchmark of k-mer counting and PORT-EK-v2 stages.
# usage:  bash scripts/s05_benchmark.sh FASTA K THREADS REPS PROJECT_DIR
#   FASTA        one concatenated FASTA of all sequences used for the counters (use the SAME sequences as PROJECT_DIR/input)
#   PROJECT_DIR  a PORT-EK project (config.yaml + input/). It is OVERWRITTEN by find_enriched/map/tree: use a copy.
#   env: STAGES="count_index find_enriched map tree_n250" (default: all)   PORTEK_DIR=/path/to/PORT-EK-version-2
# Output: bench_<time>.tsv (+ .logs/ with every command's output).  Rows with exit_code != 0 are FAILED runs (ignore their times).
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ $# -ge 5 ] || { echo "usage: bash $0 FASTA K THREADS REPS PROJECT_DIR"; exit 1; }
IN=$1; K=$2; T=$3; R=$4; PROJ=$5; STAGES=${STAGES:-"count_index find_enriched map tree_n250"}
[ -f "$IN" ] || { echo "FASTA not found: $IN"; exit 1; }
[ -f "$PROJ/config.yaml" ] || { echo "PROJECT_DIR must contain config.yaml: $PROJ"; exit 1; }
OUT=bench_$(date +%s).tsv; B="python $SD/bench_run.py $OUT"
{ echo "# $(date)"; lscpu | grep -E "Model name|^CPU\(s\)"; free -g | head -2; python -V; jellyfish --version; kmc 2>&1 | head -1; } > ${OUT%.tsv}_environment.txt 2>&1
have() { command -v "$1" >/dev/null 2>&1 || { echo "WARNING: $1 not found in PATH -> skipped"; return 1; }; }
for REP in $(seq $R); do
  if have jellyfish; then $B jellyfish count $REP -- jellyfish count -m $K -s 100M -t $T -o /tmp/jf_$$.jf $IN; fi   # strand-specific like PORT-EK (no -C); -s = hash size, which dominates peak memory
  if have kmc; then mkdir -p /tmp/kmc_$$; $B kmc count $REP -- kmc -k$K -t$T -ci1 -b -fm $IN /tmp/kmc_$$/db /tmp/kmc_$$; fi
  for S in $STAGES; do
    case $S in
      count_index)   $B portek count_index $REP -- python $SD/portek_patched_cli.py find_k $PROJ --min_k $K --max_k $K --n_jobs $T ;;
      find_enriched) $B portek find_enriched $REP -- python $SD/portek_patched_cli.py find_enriched $PROJ -k $K --n_jobs $T ;;
      map)           $B portek map $REP -- python $SD/portek_patched_cli.py map $PROJ -k $K -d 3 ;;
      tree_n250)     $B portek tree_n250 $REP -- python $SD/portek_patched_cli.py tree $PROJ -k $K -n 250 -b ;;
    esac
  done
done
python - "$OUT" <<'PY'
import sys, pandas as pd
d = pd.read_csv(sys.argv[1], sep="\t"); bad = d[d.exit_code != 0]
if len(bad): print(f"\nWARNING: {len(bad)} failed run(s) excluded from the summary (see .logs):\n", bad.groupby(["tool", "task"]).size())
g = d[d.exit_code == 0].groupby(["tool", "task"]).agg(n=("wall_s", "size"), wall_mean=("wall_s", "mean"), wall_sd=("wall_s", "std"), peak_MB_max=("peak_MB", "max")).round(2)
print("\nSummary of successful runs:\n", g if len(g) else "none"); g.to_csv(sys.argv[1].replace(".tsv", "_summary.tsv"), sep="\t")
PY
echo "raw table: $OUT"
