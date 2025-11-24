"""
count_kbins.py

Compute and print the number of linear k-bins that bskit will use,
given your kmin, kmax, and dk settings, plus the recommended ENDI index
(i.e. bin_count - 1) for use as the upper loop bound.
"""

import bskit

# ─── HYPERPARAMETERS (match your shell script) ────────────────────────────────
KMIN = 0.01        # export KMIN=0.01
KMAX = 0.4         # export KMAX=0.4
DK = 0.02          # export DK=0.02
NUM_LOW_K_BINS = 0 # default unless you’ve separately specified low-k bins
DK_HIGH = -1       # default (no separate high-k bin spacing)
# ──────────────────────────────────────────────────────────────────────────────

def main():
    edges = bskit.generate_bin_edge_list(
        kmin=KMIN,
        kmax=KMAX,
        dk=DK,
        num_lowk_bins=NUM_LOW_K_BINS,
        dk_high=DK_HIGH
    )
    bin_count = len(edges)
    print(f"Number of linear k bins: {bin_count}")
    print(f"Recommended ENDI = {bin_count - 1}")

if __name__ == "__main__":
    main()