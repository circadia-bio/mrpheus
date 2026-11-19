#!/usr/bin/env python3
"""
fetch_yasa_models.py
--------------------
Converts all eight pre-trained YASA sleep staging classifiers to gzip-compressed
LightGBM text files in inst/models/, one per channel set, with and without the
demographic features:

    yasa_eeg.txt.gz            yasa_eeg_demo.txt.gz
    yasa_eeg_eog.txt.gz        yasa_eeg_eog_demo.txt.gz
    yasa_eeg_emg.txt.gz        yasa_eeg_emg_demo.txt.gz
    yasa_eeg_eog_emg.txt.gz    yasa_eeg_eog_emg_demo.txt.gz

mrpheus::stage_epochs() picks the file that matches the channels it is given,
exactly as yasa.SleepStaging does. Each file is 3.2-3.3 MB (about 26 MB in
total). The text format is lossless: the eeg+eog+emg file written by this script
is byte-identical (md5) to the single model earlier versions of mrpheus shipped
as inst/models/yasa_staging.txt, which can be deleted.

Usage (from the package root)
-----------------------------
    pip install yasa lightgbm      # yasa >= 0.8
    python3 data-raw/fetch_yasa_models.py

Licence note
------------
YASA is released under the BSD 3-Clause License. Redistribution of the trained
model weights with attribution is consistent with this licence:
https://github.com/raphaelvallat/yasa/blob/master/LICENSE
"""
import glob
import gzip
import os
import re
import shutil
import warnings

import joblib
import yasa

warnings.filterwarnings("ignore")

OUT_DIR = os.path.join("inst", "models")


def main():
    clf_dir = os.path.join(os.path.dirname(yasa.__file__), "classifiers")
    os.makedirs(OUT_DIR, exist_ok=True)
    files = sorted(glob.glob(os.path.join(clf_dir, "clf_*_lgb_*.joblib")))
    if len(files) != 8:
        raise SystemExit(f"Expected 8 classifiers in {clf_dir}, found {len(files)}.")
    print(f"yasa {yasa.__version__}: {len(files)} classifiers")
    for path in files:
        base = os.path.basename(path)
        key = re.match(r"clf_(.+)_lgb_[\d.]+\.joblib", base).group(1)         # e.g. eeg+eog+emg+demo
        name = "yasa_" + key.replace("+", "_") + ".txt"
        txt = os.path.join(OUT_DIR, name)
        joblib.load(path).booster_.save_model(txt)
        with open(txt, "rb") as fi, gzip.open(txt + ".gz", "wb", compresslevel=9) as fo:
            shutil.copyfileobj(fi, fo)
        os.remove(txt)
        print(f"  {base:40s} -> {name}.gz ({os.path.getsize(txt + '.gz') / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
