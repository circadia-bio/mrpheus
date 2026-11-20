#!/usr/bin/env python3
"""
generate_yasa_bandpower_parity_fixture.py
-----------------------------------------
Writes tests/testthat/fixtures/yasa_bandpower_parity.csv: the output of
yasa.bandpower() with its default arguments (absolute and relative band power
and TotalAbsPow, one 30 s epoch at a time) on a deterministic synthetic EEG
channel at 100 Hz and at 256 Hz (no participant data).

tests/testthat/test-yasa-bandpower-parity.R generates the same signal from the
same formulas and requires compute_band_power() to reproduce the table.

Signal definition (micro-volts; the same EEG channel as
generate_yasa_parity_fixture.py, with the sampling rate as a parameter):
  noise[i] : Park-Miller minimal-standard LCG, state <- 16807 * state mod
             (2^31 - 1), u = state / (2^31 - 1) - 0.5, seed 12345. 16807 * state
             stays below 2^53, so the stream is exact in double precision.
  t        : i / fs
  env      : 1 + 0.6 sin(2 pi t / 1200)
  EEG      : 30 env (sin(2 pi 1.1 t + 0.3) + 0.5 sin(2 pi 10.3 t)
                     + 0.3 sin(2 pi 2.7 t + 1)) + 12 noise

Usage:  python data-raw/generate_yasa_bandpower_parity_fixture.py
Needs:  yasa >= 0.8, numpy, pandas.
"""
import os
import warnings

import numpy as np
import pandas as pd
import yasa

warnings.filterwarnings("ignore")

N_EPOCHS = 6
RATES = (100, 256)
BANDS = ["delta", "theta", "alpha", "sigma", "beta", "gamma"]


def lcg(n, seed):
    out = np.empty(n)
    state = seed
    for i in range(n):
        state = (16807 * state) % 2147483647
        out[i] = state / 2147483647 - 0.5
    return out


def eeg(fs):
    n = fs * 30 * N_EPOCHS
    t = np.arange(n) / fs
    tp = 2 * np.pi
    env = 1 + 0.6 * np.sin(tp * t / 1200)
    return 30 * env * (np.sin(tp * 1.1 * t + 0.3) + 0.5 * np.sin(tp * 10.3 * t)
                       + 0.3 * np.sin(tp * 2.7 * t + 1)) + 12 * lcg(n, 12345)


if __name__ == "__main__":
    rows = []
    for fs in RATES:
        x = eeg(fs)
        print(f"fs {fs}: first two samples {x[0]!r}, {x[1]!r}")
        spe = fs * 30
        for e in range(N_EPOCHS):
            seg = x[e * spe:(e + 1) * spe]
            absolute = yasa.bandpower(seg, fs, ch_names=["EEG"], relative=False)
            relative = yasa.bandpower(seg, fs, ch_names=["EEG"], relative=True)
            row = {"fs": fs, "epoch": e + 1}
            row.update({b: absolute[b.capitalize()].iloc[0] for b in BANDS})
            row.update({f"rel_{b}": relative[b.capitalize()].iloc[0] for b in BANDS})
            row["total_power"] = absolute["TotalAbsPow"].iloc[0]
            rows.append(row)
    out = os.path.join("tests", "testthat", "fixtures", "yasa_bandpower_parity.csv")
    pd.DataFrame(rows).to_csv(out, index=False, float_format="%.17g")
    print("written", out, f"(yasa {yasa.__version__})")
