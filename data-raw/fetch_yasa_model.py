#!/usr/bin/env python3
"""
fetch_yasa_model.py
-------------------
Extracts the pre-trained LightGBM sleep staging model from YASA and saves it
as a plain-text file for use in mrpheus.

The LightGBM model format (.txt) is fully cross-language: a model trained in
Python loads identically in R's lightgbm package (both are thin wrappers around
the same C++ library). Once saved here it is bundled in inst/models/ and ships
with the package — no Python required by end users.

Usage
-----
    python3 -m venv /tmp/yasa_env
    source /tmp/yasa_env/bin/activate
    pip install yasa lightgbm
    python3 data-raw/fetch_yasa_model.py

Requirements
------------
    yasa >= 0.6.0
    lightgbm >= 4.0.0
    joblib (installed as a yasa dependency)

Model selection
---------------
Uses clf_eeg+eog+emg_lgb_0.5.0.joblib — the most complete channel set
(EEG + EOG + EMG). This is the model called by stage_epochs() in mrpheus.

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

import pathlib

try:
    import joblib
    import yasa
except ImportError as e:
    raise ImportError(
        "Please install dependencies first:\n"
        "    pip install yasa lightgbm"
    ) from e

# Model file to extract (EEG + EOG + EMG, v0.5.0)
MODEL_NAME = "clf_eeg+eog+emg_lgb_0.5.0.joblib"

OUT_DIR  = pathlib.Path(__file__).parent.parent / "inst" / "models"
OUT_PATH = OUT_DIR / "yasa_staging.txt"

def find_booster(obj):
    """Recursively locate the underlying lgb.Booster in a sklearn estimator."""
    # LGBMClassifier exposes .booster_
    if hasattr(obj, "booster_"):
        return obj.booster_
    # Pipeline: walk steps
    if hasattr(obj, "steps"):
        for _, step in obj.steps:
            result = find_booster(step)
            if result is not None:
                return result
    # Attribute named clf, estimator, or _clf
    for attr in ("clf", "estimator", "_clf", "_estimator"):
        if hasattr(obj, attr):
            result = find_booster(getattr(obj, attr))
            if result is not None:
                return result
    return None


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    yasa_dir   = pathlib.Path(yasa.__file__).parent
    model_path = yasa_dir / "classifiers" / MODEL_NAME

    if not model_path.exists():
        # List available models to help the user pick an alternative
        available = sorted((yasa_dir / "classifiers").glob("*.joblib"))
        names     = [f.name for f in available]
        raise FileNotFoundError(
            f"Expected model not found: {model_path}\n"
            f"Available classifiers:\n  " + "\n  ".join(names)
        )

    print(f"Loading: {model_path.name} "
          f"({model_path.stat().st_size / 1024:.0f} KB)")

    clf = joblib.load(model_path)

    booster = find_booster(clf)
    if booster is None:
        raise RuntimeError(
            f"Could not locate lgb.Booster inside loaded object.\n"
            f"Object type: {type(clf)}\n"
            f"Attributes:  {[a for a in dir(clf) if not a.startswith('__')]}"
        )

    booster.save_model(str(OUT_PATH))

    size_kb = OUT_PATH.stat().st_size / 1024
    print(f"Saved:  {OUT_PATH}  ({size_kb:.0f} KB)")
    print(f"YASA version: {yasa.__version__}")
    print("\nNext: commit inst/models/yasa_staging.txt to the package repo.")


if __name__ == "__main__":
    main()
