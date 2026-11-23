#!/usr/bin/env python3
"""R1-Major6: show how k depends on L (gapped vs ACGT-only) and on the knee rule; then sensitivity over k.
Usage: python s04_k_selection_and_sensitivity.py L_gapped L_acgt          e.g. 19918 8995 (DNA)   12621 8889 (RNA)
Compute L with:  python - <<'X'  (mean len(seq) vs mean count of ACGT per FASTA record)  X
Sensitivity over k: bash s04b_k_sensitivity.sh PROJECT 13 15 17 19 21
"""
import sys, numpy as np
Lg, La = float(sys.argv[1]), float(sys.argv[2]); ks = np.arange(5, 32, 2)
def spec(L): return np.exp(-(L - ks + 1) ** 2 / (2 * 4.0 ** ks))
def chord(y, orig=False):
    x = (ks - ks[0]) / (ks[-1] - ks[0]); y = (y - y.min()) / (y.max() - y.min()); dx, dy = x[-1] - x[0], y[-1] - y[0]
    t = dy * x - dx * y + x[-1] * y[0] - x[0] * y[-1]
    if orig: t = dy * x - dx * y + x[-1] * y[0] - dx * y + x[-1] * y[0] - x[0] * y[-1]   # original (duplicated terms)
    return ks[np.argmax(np.abs(t))]
for name, L in (("gapped L (as run)", Lg), ("ACGT-only L", La)):
    y = spec(L); print(f"{name:20s} L={L:7.0f}  original rule -> k={chord(y, True)}   corrected chord -> k={chord(y)}   P(no collision) at k=13/15/17: {y[4]:.3f}/{y[5]:.3f}/{y[6]:.3f}")
