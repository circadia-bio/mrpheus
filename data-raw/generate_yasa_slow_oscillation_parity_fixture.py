#!/usr/bin/env python3
"""
generate_yasa_slow_oscillation_parity_fixture.py
------------------------------------------------
Writes fixtures from YASA on synthetic data (no participant data):

  tests/testthat/fixtures/yasa_sw_events.csv   the events of yasa.sw_detect()
        (default arguments) at 100 Hz and at 256 Hz, with a stage mask that
        excludes some epochs and cuts through events
  tests/testthat/fixtures/yasa_sw_filter.csv   taps of the FIR band-pass that
        sw_detect applies (mne.filter.create_filter, 0.3-1.5 Hz, 0.2 Hz
        transition bands) at 100 and 256 Hz, every 25th tap
  tests/testthat/fixtures/yasa_sw_crafted_signal.csv   a crafted, already
        band-passed signal (integer uV, 100 Hz) and a stage mask, built to reach
        the rules of sw_detect that a smooth signal does not: peak plateaus,
        double-humped waves (events that share a start), durations and
        amplitudes exactly on a limit, exact zeros between waves, a wave before
        the first zero crossing, positive waves too small to count as peaks
  tests/testthat/fixtures/yasa_sw_crafted_events.csv   the events YASA's own
        detection code finds in it (sw_detect with its band-pass replaced by the
        identity)

tests/testthat/test-yasa-slow-oscillation-parity.R generates the synthetic EEG
from the same formulas and requires compute_slow_oscillations() to reproduce
the events, mrpheus' filter design to reproduce the taps, and the detection to
reproduce the events of the crafted signal.

Signal definition (micro-volts; 10 epochs of 30 s; mirrored in the R test):
  noise[i] : Park-Miller minimal-standard LCG, state <- 16807 * state mod
             (2^31 - 1), u = state / (2^31 - 1) - 0.5, seed 12345
  t        : i / fs
  env      : 1 + 0.5 sin(2 pi t / 180)
  EEG      : 70 env (sin(2 pi 0.75 t) + 0.45 sin(2 pi 1.15 t + 0.5))
             + 14 sin(2 pi 0.3 t + 2) + 18 noise
Stage mask : epochs 1-3, 5-7 and 9-10 are accepted (epochs 4 and 8 are not).

Usage:  python data-raw/generate_yasa_slow_oscillation_parity_fixture.py
Needs:  yasa >= 0.8, mne, numpy, pandas.
"""
import os
import warnings

import mne
import numpy as np
import pandas as pd
import yasa

warnings.filterwarnings("ignore")
mne.set_log_level("ERROR")

N_EPOCHS = 10
RATES = (100, 256)
COLS = ["Start", "NegPeak", "MidCrossing", "PosPeak", "End", "Duration",
        "ValNegPeak", "ValPosPeak", "PTP", "Slope", "Frequency"]


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
    env = 1 + 0.5 * np.sin(tp * t / 180)
    return (70 * env * (np.sin(tp * 0.75 * t) + 0.45 * np.sin(tp * 1.15 * t + 0.5))
            + 14 * np.sin(tp * 0.3 * t + 2) + 18 * lcg(n, 12345))


def stage_mask(fs):
    flag = (np.arange(1, N_EPOCHS + 1) % 4) != 0
    return np.repeat(flag, fs * 30)


def crafted(sf=100, seed=2024, n_lobes=150):
    """A band-passed-looking signal made of alternating negative and positive
    lobes (integer uV), with a stage mask."""
    rng = np.random.default_rng(seed)
    x = []
    lobes = []                       # (sign, length, amplitude, humps, zero_before)
    sign = -1                        # starts with a negative lobe: no zero crossing before its peak
    # deliberate pairs (negative amplitude, negative length, positive amplitude,
    # positive length in samples) with a limit exactly on the boundary, right
    # after the first wave: peak-to-peak 75 / 350 (exclusive), negative duration
    # 0.3 / 1.5 s and positive duration 0.1 / 1.0 s (exclusive), each next to
    # the nearest accepted value
    deliberate = [(40, 31, 35, 51), (41, 31, 35, 51), (200, 61, 150, 61), (200, 61, 149, 61),
                  (60, 30, 60, 51), (60, 31, 60, 51), (60, 150, 60, 51), (60, 149, 60, 51),
                  (60, 51, 60, 10), (60, 51, 60, 11), (60, 51, 60, 100), (60, 51, 60, 99)]
    for k in range(n_lobes):
        if k == 1:
            for na, nl, pa, pl in deliberate:
                lobes.append((-1, nl, na, 1, False))
                lobes.append((1, pl, pa, 1, False))
            sign = -1
        if sign < 0:
            length = int(rng.choice([30, 31, 150, 151, *range(28, 156, 1)]))
            amp = int(rng.choice([39, 40, 41, 200, 201, *range(30, 225, 5)]))
        else:
            length = int(rng.choice([10, 11, 100, 101, *range(8, 106, 1)]))
            amp = int(rng.choice([9, 10, 11, 35, 36, 150, 151, *range(6, 160, 4)]))
        humps = 2 if rng.random() < 0.15 else 1
        zero_before = rng.random() < 0.25 and k > 0
        lobes.append((sign, length, amp, humps, zero_before))
        sign = -sign
    for sign, length, amp, humps, zero_before in lobes:
        if zero_before:
            x.append(0)              # an exact zero between waves
        t = (np.arange(length) + 0.5) / length
        shape = np.sin(np.pi * t) if humps == 1 else (np.sin(np.pi * t) + 0.6 * np.sin(3 * np.pi * t)) / 1.0
        shape = shape / shape.max()
        vals = np.rint(sign * amp * shape).astype(int)
        vals[vals * sign <= 0] = sign    # never an accidental zero or sign change
        x.extend(vals.tolist())
    x = np.array(x, dtype=float)
    mask = np.ones(len(x), dtype=bool)
    for start in range(1500, len(x) - 600, 2600):      # a few stretches outside the stage mask
        mask[start:start + int(rng.integers(150, 500))] = False
    return x, mask


def identity_sw_detect(x, sf, mask):
    """yasa.sw_detect with its band-pass replaced by the identity, so that
    YASA's own detection code runs on an already filtered signal."""
    from yasa import detection
    original = detection.filter_data
    detection.filter_data = lambda data, sf, *args, **kwargs: data
    try:
        res = yasa.sw_detect(x, sf, hypno=np.where(mask, 2, 0).astype(int), include=(2, 3))
    finally:
        detection.filter_data = original
    return res


if __name__ == "__main__":
    fixtures = os.path.join("tests", "testthat", "fixtures")
    events, taps = [], []
    for fs in RATES:
        x = eeg(fs)
        print(f"fs {fs}: first two samples {x[0]!r}, {x[1]!r}")
        hypno = np.where(stage_mask(fs), 2, 0).astype(int)
        ev = yasa.sw_detect(x, fs, hypno=hypno, include=(2, 3)).summary()[COLS].copy()
        ev.insert(0, "fs", fs)
        events.append(ev)
        print(f"fs {fs}: {len(ev)} events")
        h = mne.filter.create_filter(None, fs, 0.3, 1.5, l_trans_bandwidth=0.2,
                                     h_trans_bandwidth=0.2, method="fir",
                                     fir_design="firwin", fir_window="hamming",
                                     phase="zero", verbose=False)
        idx = np.arange(0, len(h), 25)
        taps.append(pd.DataFrame({"fs": fs, "n_taps": len(h), "tap": idx + 1, "h": h[idx]}))
    pd.concat(events).to_csv(os.path.join(fixtures, "yasa_sw_events.csv"), index=False,
                             float_format="%.15g")
    pd.concat(taps).to_csv(os.path.join(fixtures, "yasa_sw_filter.csv"), index=False,
                           float_format="%.17g")
    xc, mc = crafted()
    ev = identity_sw_detect(xc, 100, mc).summary()[COLS]
    print(f"crafted signal: {len(xc)} samples, {len(ev)} events")
    pd.DataFrame({"x": xc.astype(int), "mask": mc.astype(int)}).to_csv(
        os.path.join(fixtures, "yasa_sw_crafted_signal.csv"), index=False)
    ev.to_csv(os.path.join(fixtures, "yasa_sw_crafted_events.csv"), index=False,
              float_format="%.15g")
    print("written fixtures (yasa", yasa.__version__ + ", mne", mne.__version__ + ")")
