#!/usr/bin/env python3
"""
generate_yasa_spindle_parity_fixture.py
---------------------------------------
Writes fixtures from YASA on a deterministic synthetic EEG channel with spindle
bursts (no participant data):

  tests/testthat/fixtures/yasa_sp_events.csv    the events of yasa.spindles_detect()
        (default arguments) at 100 Hz and at 256 Hz, with a stage mask that
        excludes some epochs ("base"), and at 100 Hz two variants that reach rules
        the base signal does not: "loud" (strong 13.5 Hz activity in the excluded
        epoch 4, so a threshold computed over all samples instead of the accepted
        ones would be 4 times higher) and "x4" (the signal times 4, so that the
        moving-RMS threshold exceeds YASA's cap of 10 uV)
  tests/testthat/fixtures/yasa_sp_features.csv  YASA's intermediate features on
        the 100 Hz signal: the relative sigma power of the STFT (every 10th
        frame) and the moving correlation and RMS (every 500th sample, after
        linear interpolation)

tests/testthat/test-yasa-spindle-parity.R generates the same signal from the
same formulas and requires compute_spindles() to reproduce both.

Signal definition (micro-volts; 12 epochs of 30 s; mirrored in the R test):
  noise[i] : Park-Miller minimal-standard LCG, state <- 16807 * state mod
             (2^31 - 1), u = state / (2^31 - 1) - 0.5, seed 12345
  t        : i / fs
  background : 28 (1 + 0.4 sin(2 pi t / 150)) sin(2 pi 0.8 t + 0.4)
               + 10 sin(2 pi 6.1 t + 1) + 14 noise
  bursts   : for each (c, w, a, f) in BURSTS (see below):
               a exp(-((t - c) / w)^2) sin(2 pi f (t - c))
Stage mask : epochs 1-3, 5-7 and 9-11 are accepted (epochs 4, 8 and 12 are not).
Variants (100 Hz):  loud = signal + 120 sin(2 pi 13.5 t) (1 + 0.5 sin(2 pi 0.3 t))
                    for t in epoch 4 (90-120 s);  x4 = 4 * signal.

Usage:  python data-raw/generate_yasa_spindle_parity_fixture.py
Needs:  yasa >= 0.8, mne, numpy, pandas.
"""
import os
import warnings

import mne
import numpy as np
import pandas as pd
import yasa
from mne.filter import filter_data

warnings.filterwarnings("ignore")
mne.set_log_level("ERROR")

N_EPOCHS = 12
RATES = (100, 256)
COLS = ["Start", "Peak", "End", "Duration", "Amplitude", "AmpFiltered", "RMS",
        "AbsPower", "RelPower", "Frequency", "Oscillations", "Symmetry"]
# (centre s, half-width s, amplitude uV, frequency Hz)
BURSTS = [(12.0, 0.45, 38, 13.5), (31.0, 0.60, 34, 12.8), (47.5, 0.35, 42, 14.2), (47.5 + 0.9, 0.30, 36, 13.0),
          (75.0, 0.55, 40, 13.6), (101.0, 0.70, 33, 13.2), (101.0 + 1.6, 0.40, 35, 14.0), (128.0, 0.50, 45, 13.8),
          (151.0, 0.25, 40, 13.5), (166.0, 0.80, 37, 12.6), (196.0, 0.55, 41, 13.9), (222.0, 0.45, 39, 13.3),
          (241.0, 0.65, 36, 14.4), (263.0, 0.50, 44, 13.1), (283.0, 0.60, 38, 13.7), (305.0, 0.40, 43, 13.4),
          (329.0, 0.55, 35, 12.9), (347.0, 0.70, 40, 13.6)]


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
    x = (28 * (1 + 0.4 * np.sin(tp * t / 150)) * np.sin(tp * 0.8 * t + 0.4)
         + 10 * np.sin(tp * 6.1 * t + 1.0) + 14 * lcg(n, 12345))
    for c, w, a, f in BURSTS:
        x += a * np.exp(-(((t - c) / w) ** 2)) * np.sin(tp * f * (t - c))
    return x


def stage_mask(fs):
    flag = (np.arange(1, N_EPOCHS + 1) % 4) != 0
    return np.repeat(flag, fs * 30)


if __name__ == "__main__":
    fixtures = os.path.join("tests", "testthat", "fixtures")
    events = []
    variants = [("base", fs) for fs in RATES] + [("loud", 100), ("x4", 100)]
    for variant, fs in variants:
        x = eeg(fs)
        if variant == "loud":
            t = np.arange(x.size) / fs
            seg = slice(3 * fs * 30, 4 * fs * 30)
            x[seg] += 120 * np.sin(2 * np.pi * 13.5 * t[seg]) * (1 + 0.5 * np.sin(2 * np.pi * 0.3 * t[seg]))
        elif variant == "x4":
            x = 4 * x
        if variant == "base":
            print(f"fs {fs}: first two samples {x[0]!r}, {x[1]!r}")
        hypno = np.where(stage_mask(fs), 2, 0).astype(int)
        ev = yasa.spindles_detect(x, fs, hypno=hypno).summary()[COLS].copy()
        ev.insert(0, "variant", variant)
        ev.insert(1, "fs", fs)
        events.append(ev)
        print(f"{variant} fs {fs}: {len(ev)} events")
    pd.concat(events).to_csv(os.path.join(fixtures, "yasa_sp_events.csv"), index=False,
                             float_format="%.15g")

    # intermediate features at 100 Hz, from YASA's own functions
    fs = 100
    x = eeg(fs)
    xb = filter_data(x[None, :], fs, 1, 30, method="fir", verbose=False)[0]
    xs = filter_data(x[None, :], fs, 12, 15, l_trans_bandwidth=1.5, h_trans_bandwidth=1.5,
                     method="fir", verbose=False)[0]
    f, t, Sxx = yasa.stft_power(xb, fs, window=2, step=0.2, band=(1, 30), interp=False, norm=False)
    idx_sigma = np.logical_and(f >= 12, f <= 15)
    rel = Sxx[idx_sigma].sum(0) / Sxx.sum(0)
    _, mcorr = yasa.moving_transform(x=xs, y=xb, sf=fs, window=0.3, step=0.1, method="corr", interp=True)
    _, mrms = yasa.moving_transform(x=xs, sf=fs, window=0.3, step=0.1, method="rms", interp=True)
    fr = np.arange(0, t.size, 10)                       # every 10th STFT frame
    frames = pd.DataFrame({"kind": "rel_pow", "index": fr + 1, "time": t[fr], "value": rel[fr]})
    sel = np.arange(0, x.size, 500)
    samples = pd.concat([
        pd.DataFrame({"kind": "corr", "index": sel + 1, "time": sel / fs, "value": mcorr[sel]}),
        pd.DataFrame({"kind": "rms", "index": sel + 1, "time": sel / fs, "value": mrms[sel]}),
    ])
    pd.concat([frames, samples]).to_csv(os.path.join(fixtures, "yasa_sp_features.csv"), index=False,
                                        float_format="%.15g")
    print("written fixtures (yasa", yasa.__version__ + ")")
