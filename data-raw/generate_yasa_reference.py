#!/usr/bin/env python3
"""
generate_yasa_reference.py
--------------------------
Runs YASA's internal feature extraction on a PSG EDF recording and saves the
raw 149-feature matrix as CSV. Used by validate_feature_parity.R to check that
mrpheus reproduces the same features as the Python reference implementation.

Usage
-----
    source /tmp/yasa_env/bin/activate      # same env as fetch_yasa_model.py
    python3 data-raw/generate_yasa_reference.py \\
        --edf path/to/recording.edf \\
        --eeg "EEG Fpz-Cz" \\
        --eog "EOG horizontal" \\
        --emg "EMG submental" \\
        --out data-raw/yasa_reference_features.csv

    # EOG and EMG are optional; omit or leave blank to run EEG-only.

Requirements
------------
    mne >= 0.24
    yasa >= 0.6.0
    pandas >= 1.3
"""

import argparse
import pathlib
import warnings

warnings.filterwarnings("ignore")

try:
    import mne
    import yasa
    import pandas as pd
except ImportError as e:
    raise ImportError(
        "Please install dependencies first:\n"
        "    pip install mne yasa pandas"
    ) from e


def main():
    parser = argparse.ArgumentParser(
        description="Export YASA feature matrix to CSV for mrpheus parity testing"
    )
    parser.add_argument("--edf", required=True,
                        help="Path to PSG EDF recording")
    parser.add_argument("--eeg", required=True,
                        help="EEG channel label (as it appears in the EDF header)")
    parser.add_argument("--eog", default="",
                        help="EOG channel label (optional)")
    parser.add_argument("--emg", default="",
                        help="EMG channel label (optional)")
    parser.add_argument("--out",
                        default=str(
                            pathlib.Path(__file__).parent
                            / "yasa_reference_features.csv"
                        ),
                        help="Output CSV path")
    args = parser.parse_args()

    # ── Load EDF ──────────────────────────────────────────────────────────────
    print(f"Loading: {args.edf}")
    raw = mne.io.read_raw_edf(args.edf, preload=True, verbose=False)
    raw.del_proj()
    print(f"  Duration : {raw.times[-1] / 3600:.2f} h")
    print(f"  Channels : {raw.ch_names}")

    # ── Initialise SleepStaging ───────────────────────────────────────────────
    kwargs = {"eeg_name": args.eeg}
    if args.eog:
        kwargs["eog_name"] = args.eog
    if args.emg:
        kwargs["emg_name"] = args.emg

    print(f"\nSleepStaging channels:")
    print(f"  EEG = {args.eeg!r}")
    print(f"  EOG = {args.eog!r or '(none)'}")
    print(f"  EMG = {args.emg!r or '(none)'}")

    sls = yasa.SleepStaging(raw, **kwargs)

    # ── Extract features ──────────────────────────────────────────────────────
    # predict() triggers _compute_features() and stores the matrix in
    # sls._features (pandas DataFrame: n_epochs × 149).
    print("\nExtracting features (this takes a moment)...")
    _ = sls.predict()

    features: pd.DataFrame = sls._features
    print(f"\nFeature matrix : {features.shape[0]} epochs × {features.shape[1]} features")
    print(f"First 5 column names : {list(features.columns[:5])}")
    print(f"Last  5 column names : {list(features.columns[-5:])}")

    # Basic sanity: should be 149 columns and match the model feature_names
    if features.shape[1] != 149:
        print(f"\nWARNING: expected 149 features, got {features.shape[1]}")

    # ── Save ──────────────────────────────────────────────────────────────────
    out_path = pathlib.Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(out_path, index=False)

    print(f"\nSaved  : {out_path}")
    print(f"Size   : {out_path.stat().st_size / 1024:.1f} KB")
    print(f"\nYASA version : {yasa.__version__}")
    print(f"MNE  version : {mne.__version__}")
    print("\nNext step: run data-raw/validate_feature_parity.R")


if __name__ == "__main__":
    main()
