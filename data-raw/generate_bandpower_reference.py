#!/usr/bin/env python3
"""
generate_bandpower_reference.py
--------------------------------
Runs YASA's bandpower() on a PSG EDF recording epoch-by-epoch and saves the
results as CSV. Used by validate_band_power_parity.R to check that
compute_band_power() in mrpheus reproduces the same values as the Python
reference.

The Welch parameters match the notebook defaults:
    win_sec=4, noverlap=200, nfft=1024, window='hann',
    detrend='constant', average='mean'

Usage
-----
    source /tmp/yasa_env/bin/activate

    # With the YASA/PhysioNet example file (downloads automatically):
    python3 data-raw/generate_bandpower_reference.py --download-example

    # With your own EDF:
    python3 data-raw/generate_bandpower_reference.py \\
        --edf path/to/recording.edf \\
        --channel "EEG Fpz-Cz" \\
        --sf 100 \\
        --out data-raw/yasa_reference_bandpower.csv

Requirements
------------
    mne >= 0.24, yasa >= 0.6.0, numpy, pandas, scipy
"""

import argparse
import pathlib
import warnings

warnings.filterwarnings("ignore")

try:
    import mne
    import numpy as np
    import pandas as pd
    import yasa
except ImportError as e:
    raise ImportError(
        "Please install dependencies:\n"
        "    pip install mne yasa pandas numpy"
    ) from e

# ── Band definitions matching compute_band_power() defaults ──────────────────
BANDS = [
    (0.5,  4,  "delta"),
    (4,    8,  "theta"),
    (8,   12,  "alpha"),
    (12,  16,  "sigma"),
    (16,  30,  "beta"),
    (30,  45,  "gamma"),
]

WELCH_KWARGS = {
    "noverlap": 200,
    "nfft":     1024,
    "window":   "hann",
    "average":  "mean",
    "detrend":  "constant",
}

EPOCH_S    = 30
WIN_SEC    = 4


def fetch_physionet_example():
    from mne.datasets.sleep_physionet.age import fetch_data
    paths = fetch_data(subjects=[0], recording=[1], on_missing="warn")
    return pathlib.Path(paths[0][0])


def bandpower_per_epoch(raw, channel, epoch_s=30, relative=True):
    """
    Segment a raw MNE recording into epochs and compute YASA bandpower for each.
    Returns a DataFrame with columns: epoch, channel, band, power, relative_power.
    """
    sf       = raw.info["sfreq"]
    data, _  = raw[channel]                  # (1, n_samples), MNE returns Volts
    signal   = data[0] * 1e6                 # convert V -> uV to match edfReader
    n_samp   = int(epoch_s * sf)
    n_epochs = int(len(signal) // n_samp)

    rows = []
    for i in range(n_epochs):
        seg = signal[i * n_samp : (i + 1) * n_samp]

        # YASA bandpower expects (n_channels, n_times) or (n_times,)
        bp_abs = yasa.bandpower(
            seg,
            sf          = sf,
            win_sec     = WIN_SEC,
            relative    = False,
            bandpass    = False,
            bands       = BANDS,
            kwargs_welch = WELCH_KWARGS,
        )
        bp_rel = yasa.bandpower(
            seg,
            sf          = sf,
            win_sec     = WIN_SEC,
            relative    = True,
            bandpass    = False,
            bands       = BANDS,
            kwargs_welch = WELCH_KWARGS,
        )

        for fmin, fmax, band_name in BANDS:
            col = band_name  # YASA preserves the name as-is
            rows.append({
                "epoch":          i + 1,          # 1-indexed to match R
                "channel":        channel,
                "band":           band_name,
                "power":          float(bp_abs[col].iloc[0]),
                "relative_power": float(bp_rel[col].iloc[0]),
                "total_power":    float(bp_abs[[b[2] for b in BANDS]
                                         ].sum(axis=1).iloc[0]),
            })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edf",     default=None, help="Path to EDF file")
    parser.add_argument("--channel", default="EEG Fpz-Cz",
                        help="Channel to analyse (default: 'EEG Fpz-Cz')")
    parser.add_argument("--out",     default="data-raw/yasa_reference_bandpower.csv",
                        help="Output CSV path")
    parser.add_argument("--download-example", action="store_true",
                        help="Download SC4001E0-PSG.edf from PhysioNet via MNE")
    parser.add_argument("--n-epochs", type=int, default=None,
                        help="Limit to first N epochs (default: all)")
    args = parser.parse_args()

    if args.download_example:
        edf_path = fetch_physionet_example()
        print(f"Using example EDF: {edf_path}")
    elif args.edf:
        edf_path = pathlib.Path(args.edf)
    else:
        parser.error("Supply --edf or --download-example")

    print(f"Loading: {edf_path}")
    raw = mne.io.read_raw_edf(edf_path, preload=True, verbose=False)
    raw.resample(100, verbose=False)     # match mrpheus default staging sr

    if args.channel not in raw.ch_names:
        raise ValueError(
            f"Channel '{args.channel}' not found.\n"
            f"Available: {raw.ch_names}"
        )

    print(f"Channel : {args.channel}  |  sf = {raw.info['sfreq']} Hz")

    df = bandpower_per_epoch(raw, args.channel, epoch_s=EPOCH_S)
    if args.n_epochs:
        df = df[df["epoch"] <= args.n_epochs]

    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"Saved {len(df)} rows → {out_path}")
    print(df.groupby("band")[["power", "relative_power"]].mean().round(6))


if __name__ == "__main__":
    main()
