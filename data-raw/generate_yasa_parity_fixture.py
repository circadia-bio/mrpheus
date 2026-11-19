#!/usr/bin/env python3
"""
generate_yasa_parity_fixture.py
-------------------------------
Writes two fixtures from YASA on a deterministic synthetic three-channel
recording (no participant data):

  tests/testthat/fixtures/yasa_parity_features.csv       the 149 features YASA
        computes for EEG + EOG + EMG (float32 precision), the primary guard
  tests/testthat/fixtures/yasa_parity_probabilities.csv  posterior probabilities
        for all eight classifier configurations, an end-to-end check

tests/testthat/test-yasa-parity.R generates the same recording from the same
formulas and requires mrpheus to reproduce both.

Signal definition (100 Hz, 30 epochs = 15 min; mirrored in the R test):
  noise_k[i]  : Park-Miller minimal-standard LCG, state <- 16807 * state mod
                (2^31 - 1), u = state / (2^31 - 1) - 0.5, one stream per channel
                (seeds 12345, 67890, 424242). 16807 * state stays below 2^53, so
                the stream is exact in double precision in both languages.
  t           : i / 100
  env         : 1 + 0.6 sin(2 pi t / 1200)
  EEG         : 30 env (sin(2 pi 1.1 t + 0.3) + 0.5 sin(2 pi 10.3 t)
                        + 0.3 sin(2 pi 2.7 t + 1)) + 12 noise_eeg
  EOG         : 20 (sin(2 pi 0.4 t) (1 + 0.8 sin(2 pi t / 900))
                    + 0.3 sin(2 pi 5.1 t)) + 8 noise_eog
  EMG         : 4 (1 + 0.9 sin(2 pi t / 700 + 2)) 2 noise_emg + 2 sin(2 pi 13.7 t)
  (micro-volts)

Classifiers: eeg, eeg+eog, eeg+emg, eeg+eog+emg, each with and without
metadata (age 50, male 1).

Usage:  python data-raw/generate_yasa_parity_fixture.py
Needs:  yasa >= 0.8, mne, numpy, pandas.
"""
import itertools
import os
import warnings

import mne
import numpy as np
import pandas as pd
import yasa

warnings.filterwarnings("ignore")
mne.set_log_level("ERROR")

FS, N_EPOCHS = 100, 30
AGE, MALE = 50, 1


def lcg(n, seed):
    out = np.empty(n)
    state = seed
    for i in range(n):
        state = (16807 * state) % 2147483647
        out[i] = state / 2147483647 - 0.5
    return out


def signals():
    n = FS * 30 * N_EPOCHS
    t = np.arange(n) / FS
    tp = 2 * np.pi
    env = 1 + 0.6 * np.sin(tp * t / 1200)
    eeg = 30 * env * (np.sin(tp * 1.1 * t + 0.3) + 0.5 * np.sin(tp * 10.3 * t)
                      + 0.3 * np.sin(tp * 2.7 * t + 1)) + 12 * lcg(n, 12345)
    eog = 20 * (np.sin(tp * 0.4 * t) * (1 + 0.8 * np.sin(tp * t / 900))
                + 0.3 * np.sin(tp * 5.1 * t)) + 8 * lcg(n, 67890)
    emg = 4 * (1 + 0.9 * np.sin(tp * t / 700 + 2)) * 2 * lcg(n, 424242) + 2 * np.sin(tp * 13.7 * t)
    return eeg, eog, emg


if __name__ == "__main__":
    eeg, eog, emg = signals()
    info = mne.create_info(["EEG", "EOG", "EMG"], FS, ["eeg", "eog", "emg"])
    raw = mne.io.RawArray(np.vstack([eeg, eog, emg]) * 1e-6, info, verbose=False)   # volts

    rows = []
    for use_eog, use_emg, demo in itertools.product((False, True), (False, True), (False, True)):
        key = "eeg" + ("_eog" if use_eog else "") + ("_emg" if use_emg else "") + ("_demo" if demo else "")
        sls = yasa.SleepStaging(raw, eeg_name="EEG",
                                eog_name="EOG" if use_eog else None,
                                emg_name="EMG" if use_emg else None,
                                metadata={"age": AGE, "male": MALE} if demo else None)
        proba = sls.predict().proba.copy()
        proba.columns = ["W", "N1", "N2", "N3", "REM"]
        for epoch, row in enumerate(proba[["N1", "N2", "N3", "REM", "W"]].itertuples(index=False), start=1):
            rows.append((key, epoch, *row))

    feats = yasa.SleepStaging(raw, eeg_name="EEG", eog_name="EOG", emg_name="EMG").get_features()
    fpath = os.path.join("tests", "testthat", "fixtures", "yasa_parity_features.csv")
    feats.to_csv(fpath, float_format="%.9g")                       # index is the 0-based epoch
    print(f"written {fpath} ({feats.shape[0]} epochs x {feats.shape[1]} features)")

    out = pd.DataFrame(rows, columns=["classifier", "epoch", "prob_N1", "prob_N2", "prob_N3", "prob_REM", "prob_W"])
    path = os.path.join("tests", "testthat", "fixtures", "yasa_parity_probabilities.csv")
    out.to_csv(path, index=False, float_format="%.12g")
    print(f"written {path} ({len(out)} rows; yasa {yasa.__version__})")
    print("first samples (generator check):", eeg[:2].tolist(), eog[:2].tolist(), emg[:2].tolist())
