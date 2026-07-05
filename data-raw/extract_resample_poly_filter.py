#!/usr/bin/env python3
"""
extract_resample_poly_filter.py
--------------------------------
Extracts the internal FIR anti-aliasing filter that scipy.signal.resample_poly
uses for upsampling from 1 Hz to 100 Hz, and saves it for use in mrpheus.

scipy.signal.resample_poly(x, up=100, down=1) internally designs a
Kaiser-windowed lowpass FIR via:

    half_len = n_zeros * max(up, down)       # default n_zeros = 10
    n        = 2 * half_len + 1              # 2001 taps
    window   = ('kaiser', 5.0)               # default Kaiser beta
    b        = firwin(n, 1.0/max(up,down), window=window)  # cutoff = 0.01
    b       *= up                            # scale by up-factor

The filter is then applied via polyphase decomposition: b is split into
`up` (=100) polyphase branches, each of length ceil(n/up) = 21 taps.

This script replicates that design exactly, then saves:
  inst/filters/resample_poly_100hz.csv   — the 2001 scaled FIR coefficients
  inst/filters/resample_poly_100hz.json  — metadata (up, down, n_taps, etc.)

Usage
-----
    source /tmp/yasa_env/bin/activate
    python3 data-raw/extract_resample_poly_filter.py
"""

import json
import pathlib

try:
    import numpy as np
    import scipy.signal as sig
    import scipy
except ImportError as e:
    raise ImportError("pip install numpy scipy") from e

OUT_DIR = pathlib.Path(__file__).parent.parent / "inst" / "filters"

# ── Parameters (must match scipy.signal.resample_poly defaults) ────────────────
UP        = 100
DOWN      = 1
N_ZEROS   = 10          # scipy default: upfirdn half_len = n_zeros * max(up,down)
WINDOW    = 'boxcar'    # MNE's raw.resample() uses boxcar (not scipy's kaiser default)


def design_resample_poly_filter(up, down, n_zeros, window):
    """Replicate scipy.signal.resample_poly's internal FIR design."""
    # Simplified up/down ratio
    from math import gcd
    g   = gcd(up, down)
    up_ = up   // g
    dn_ = down // g

    half_len = n_zeros * max(up_, dn_)
    n        = 2 * half_len + 1
    cutoff   = 1.0 / max(up_, dn_)          # normalised to Nyquist of the higher rate

    b = sig.firwin(n, cutoff, window=window)
    b *= up_                                  # scale so output has correct amplitude
    return b, n, half_len, cutoff, up_, dn_


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"NumPy   : {np.__version__}")
    print(f"SciPy   : {scipy.__version__}")
    print(f"Up/down : {UP}/{DOWN}")

    b, n, half_len, cutoff, up_r, dn_r = design_resample_poly_filter(
        UP, DOWN, N_ZEROS, WINDOW
    )

    print(f"Filter length : {n} taps")
    print(f"half_len      : {half_len}")
    print(f"Cutoff (norm) : {cutoff:.6f}")
    print(f"First 5 coefs : {b[:5]}")
    print(f"Last  5 coefs : {b[-5:]}")
    print(f"Sum of coefs  : {b.sum():.6f}  (should be ≈ 1.0 before up-scaling, "
          f"≈ {up_r:.1f} after)")

    # ── Verify by comparing with scipy's own resampling ───────────────────────
    # Correct approach: use upfirdn which is what resample_poly calls internally,
    # then trim the same way (skip first n_b//2 samples = the filter half-delay).
    from scipy.signal import upfirdn
    rng  = np.random.default_rng(42)
    x    = rng.standard_normal(200)
    ref  = sig.resample_poly(x, UP, DOWN, window=WINDOW)
    full = upfirdn(b, x, up=UP, down=DOWN)   # length = N*UP + n_b - 1
    n_pre = b.size // 2                       # = half_len = 1000
    our  = full[n_pre : n_pre + len(x) * UP]
    max_diff = np.abs(ref - our).max()
    print(f"\nMax diff vs scipy reference : {max_diff:.3e}  (should be < 1e-10)")

    # ── Save ──────────────────────────────────────────────────────────────────
    coef_path = OUT_DIR / "resample_poly_100hz.csv"
    np.savetxt(coef_path, b, delimiter=",", fmt="%.18e")
    print(f"\nSaved coefficients : {coef_path}  ({n} values, "
          f"{coef_path.stat().st_size / 1024:.1f} KB)")

    meta = {
        "up":            UP,
        "down":          DOWN,
        "up_reduced":    int(up_r),
        "down_reduced":  int(dn_r),
        "n_zeros":       N_ZEROS,
        "window":        str(WINDOW),
        "n_taps":        int(n),
        "half_len":      int(half_len),
        "cutoff_norm":   float(cutoff),
        "scipy_version": scipy.__version__,
        "description": (
            "Internal FIR anti-aliasing filter from scipy.signal.resample_poly "
            f"for up={UP}, down={DOWN}. window='{WINDOW}' matching MNE raw.resample(), "
            f"{n} taps. Coefficients are already scaled by up={UP}. "
            "Apply via polyphase decomposition: split into 100 branches of "
            f"ceil({n}/100)=21 taps each."
        ),
    }
    meta_path = OUT_DIR / "resample_poly_100hz.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Saved metadata     : {meta_path}")
    print("\nNext: commit inst/filters/resample_poly_100hz.* and run "
          "devtools::load_all() in R.")


if __name__ == "__main__":
    main()
