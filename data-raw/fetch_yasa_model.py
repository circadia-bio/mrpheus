#!/usr/bin/env python3
"""
fetch_yasa_model.py
-------------------
Extracts the pre-trained LightGBM sleep staging model from YASA and saves it
as a plain-text file for use in mrpheus.

The LightGBM model format (.txt) is fully cross-language: a model trained in
Python loads identically in R's lightgbm package (both are thin wrappers around
the same C++ library). Once saved here, copy to inst/models/yasa_staging.txt.

Usage
-----
    pip install yasa lightgbm
    python data-raw/fetch_yasa_model.py

Requirements
------------
    yasa >= 0.6.0
    lightgbm >= 4.0.0

References
----------
Vallat, R., & Walker, M. P. (2021). An open-source, high-performance tool for
automated sleep staging. eLife, 10, e70092. https://doi.org/10.7554/eLife.70092

License note
------------
YASA is released under the BSD 3-Clause License. Redistribution of the trained
model weights with attribution is consistent with this license. See:
https://github.com/raphaelvallat/yasa/blob/master/LICENSE
"""

import os
import pathlib

try:
    import yasa
    import lightgbm as lgb
except ImportError as e:
    raise ImportError(
        "Please install dependencies first:\n    pip install yasa lightgbm"
    ) from e

OUT_DIR = pathlib.Path(__file__).parent.parent / "inst" / "models"
OUT_PATH = OUT_DIR / "yasa_staging.txt"

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Instantiate SleepStaging to trigger model load.
    # The model is bundled inside the yasa package under
    # yasa/classifiers/clf_eeg+eog+emg_lgb_0.5.0.joblib (version dependent).
    # We access it via the clf attribute after a dummy instantiation.
    print("Loading YASA staging model...")
    clf = yasa.SleepStaging.__new__(yasa.SleepStaging)
    clf._load_model("eeg+eog+emg")  # loads clf.clf (LGBClassifier)

    booster = clf.clf.booster_  # underlying lgb.Booster
    booster.save_model(str(OUT_PATH))

    print(f"Model saved to: {OUT_PATH}")
    print(f"File size: {OUT_PATH.stat().st_size / 1024:.1f} KB")
    print("\nNext step: add inst/models/yasa_staging.txt to your package.")
    print("Check that the file is NOT in .gitignore (it should be committed).")


if __name__ == "__main__":
    main()
