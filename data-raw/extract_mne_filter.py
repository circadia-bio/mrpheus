#!/usr/bin/env python3
"""
extract_mne_filter.py
---------------------
Extracts the exact FIR filter coefficients that MNE/YASA applies to PSG data
and saves them as a plain CSV for use in mrpheus.

YASA always resamples data to 100 Hz before filtering, then calls:
    mne.filter.filter_data(data, sfreq=100, l_freq=0.4, h_freq=30)

This script calls mne.filter.create_filter with identical parameters to
extract the coefficients without applying them. The resulting CSV is bundled
in inst/filters/ and loaded by .bandpass_filter() in R, achieving exact
filter parity with zero Python dependency at runtime.

Usage
-----
    source /tmp/yasa_env/bin/activate
    python3 data-raw/extract_mne_filter.py

Output
------
    inst/filters/mne_bandpass_100hz.csv   — filter coefficients (one per line)
    inst/filters/mne_bandpass_100hz.json  — metadata (sfreq, l_freq, h_freq,
                                            filter_length, phase, window)
"""

import json
import pathlib

try:
    import mne
    import numpy as np
except ImportError as e:
    raise ImportError(
        "Please install dependencies first:\n"
        "    pip install mne numpy"
    ) from e

# ── Filter parameters (must match YASA's SleepStaging.fit exactly) ────────────
SFREQ   = 100.0   # YASA resamples to 100 Hz before filtering
L_FREQ  = 0.4
H_FREQ  = 30.0
METHOD  = "fir"
PHASE   = "zero"
WINDOW  = "hamming"
DESIGN  = "firwin"

OUT_DIR = pathlib.Path(__file__).parent.parent / "inst" / "filters"


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"MNE version : {mne.__version__}")
    print(f"Parameters  : sfreq={SFREQ}, l_freq={L_FREQ}, h_freq={H_FREQ}, "
          f"method={METHOD!r}, phase={PHASE!r}, window={WINDOW!r}")

    # create_filter with data=None returns the coefficients without applying them
    b = mne.filter.create_filter(
        data       = None,
        sfreq      = SFREQ,
        l_freq     = L_FREQ,
        h_freq     = H_FREQ,
        method     = METHOD,
        phase      = PHASE,
        fir_window = WINDOW,
        fir_design = DESIGN,
        verbose    = False,
    )

    print(f"Filter length: {len(b)} taps")
    print(f"First 5 coefs: {b[:5]}")
    print(f"Last  5 coefs: {b[-5:]}")

    # ── Save coefficients ──────────────────────────────────────────────────────
    coef_path = OUT_DIR / "mne_bandpass_100hz.csv"
    np.savetxt(coef_path, b, delimiter=",", fmt="%.18e")
    print(f"\nSaved coefficients: {coef_path}  ({len(b)} values, "
          f"{coef_path.stat().st_size / 1024:.1f} KB)")

    # ── Save metadata ──────────────────────────────────────────────────────────
    meta = {
        "sfreq":         SFREQ,
        "l_freq":        L_FREQ,
        "h_freq":        H_FREQ,
        "filter_length": len(b),
        "method":        METHOD,
        "phase":         PHASE,
        "fir_window":    WINDOW,
        "fir_design":    DESIGN,
        "mne_version":   mne.__version__,
        "description":   (
            "Zero-phase Hamming-windowed FIR bandpass filter designed by MNE. "
            "Applied by YASA's SleepStaging.fit() before feature extraction. "
            "Load with scan/read.csv in R, apply with gsignal::filtfilt."
        ),
    }
    meta_path = OUT_DIR / "mne_bandpass_100hz.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Saved metadata    : {meta_path}")

    print("\nNext: commit inst/filters/ and run devtools::load_all() in R.")


if __name__ == "__main__":
    main()
