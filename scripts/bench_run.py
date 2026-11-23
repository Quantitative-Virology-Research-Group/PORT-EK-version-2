#!/usr/bin/env python3
"""Run ONE command and append a row to a TSV: tool, task, rep, wall_s, peak_MB, exit_code.
usage: python bench_run.py OUT.tsv TOOL TASK REP -- command arg1 arg2 ...
Peak memory = maximum resident set size of the command and its children (resource.getrusage, Linux/macOS).
A non-zero exit code means the command FAILED; its time must not be used (last stderr lines are printed)."""
import sys, time, resource, subprocess, os
out, tool, task, rep = sys.argv[1:5]; cmd = sys.argv[sys.argv.index("--") + 1:]
logdir = out + ".logs"; os.makedirs(logdir, exist_ok=True); log = f"{logdir}/{tool}_{task}_rep{rep}.log"
t0 = time.perf_counter()
with open(log, "w") as f:
    try: rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode
    except FileNotFoundError: rc = 127; f.write(f"command not found: {cmd[0]}\n")
wall = time.perf_counter() - t0; ru = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
peak = ru / 1024 if sys.platform.startswith("linux") else ru / 1048576          # Linux: KB, macOS: bytes -> MB
if not os.path.exists(out): open(out, "w").write("tool\ttask\trep\twall_s\tpeak_MB\texit_code\n")
open(out, "a").write(f"{tool}\t{task}\t{rep}\t{wall:.2f}\t{peak:.0f}\t{rc}\n")
if rc != 0: print(f"!! FAILED ({rc}): {tool} {task} rep {rep}; see {log}\n" + "".join(open(log).readlines()[-5:]), file=sys.stderr)
