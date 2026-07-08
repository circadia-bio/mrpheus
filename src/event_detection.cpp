// src/event_detection.cpp
//
// C++ hot paths for event detection algorithms:
//   roll_rms_cpp             — centered rolling RMS for spindle detection
//   detect_so_candidates_cpp — zero-crossing scan + amplitude/duration
//                              filtering for slow oscillation detection
//
// Build: Rcpp::compileAttributes()  then  devtools::load_all()

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// ── roll_rms_cpp ──────────────────────────────────────────────────────────────
// Centered rolling RMS of signal x with window width k.
// Matches zoo::rollapply(x^2, k, function(w) sqrt(mean(w)), fill=NA, align="center"):
//   positions within half a window of each edge are filled with NA.
//
// Window alignment (matches zoo centre alignment for both odd and even k):
//   left  = floor((k-1)/2)
//   right = k - 1 - left  = ceil((k-1)/2)
//   position i covers [i - left, i + right]
//
// Note: takes the raw signal (not pre-squared); squaring is done internally.
//
// [[Rcpp::export]]
NumericVector roll_rms_cpp(NumericVector x, int k) {
    int n     = x.size();
    int left  = (k - 1) / 2;
    int right = k - 1 - left;
    NumericVector out(n, NA_REAL);

    for (int i = left; i < n - right; i++) {
        double sum = 0.0;
        for (int j = i - left; j <= i + right; j++) sum += x[j] * x[j];
        out[i] = std::sqrt(sum / (double)k);
    }
    return out;
}

// ── detect_so_candidates_cpp ─────────────────────────────────────────────────
// Scan a bandpass-filtered EEG signal for slow oscillation candidates.
// Replaces the R zero-crossing + lapply loop in compute_slow_oscillations().
//
// Algorithm (mirrors Molle et al. 2002 / YASA):
//   1. Find all zero crossings (sign changes between consecutive samples).
//   2. Walk through triplets of consecutive ZCs with step 2 (overlapping
//      negative + positive half-wave pairs).
//   3. Apply duration and peak-to-peak amplitude criteria.
//   4. Return a matrix of passing candidates.
//
// Output columns (1-indexed, matching R convention):
//   [,1] neg_start — start of negative half-wave
//   [,2] neg_end   — end of negative half-wave / start of positive half-wave
//   [,3] pos_end   — end of positive half-wave
//   [,4] neg_peak  — minimum value of sig in [neg_start, neg_end]
//   [,5] pos_peak  — maximum value of sig in [neg_end, pos_end]
//   [,6] ptp       — peak-to-peak amplitude (pos_peak - neg_peak)
//
// Indices are 1-based so the caller can use them directly in R subsetting
// and time calculations (start_s = neg_start / sr, etc.).
//
// [[Rcpp::export]]
NumericMatrix detect_so_candidates_cpp(NumericVector sig,
                                        double sr,
                                        double dur_neg_min,
                                        double dur_neg_max,
                                        double dur_pos_min,
                                        double dur_pos_max,
                                        double amp_min,
                                        double amp_max) {
    int N = sig.size();

    // Step 1: zero crossings (0-indexed positions where sign changes)
    std::vector<int> zc;
    zc.reserve(N / 10);
    for (int i = 0; i < N - 1; i++) {
        int s1 = (sig[i]     > 0.0) - (sig[i]     < 0.0);
        int s2 = (sig[i + 1] > 0.0) - (sig[i + 1] < 0.0);
        if (s1 != s2) zc.push_back(i);
    }

    int nzc = (int)zc.size();
    if (nzc < 4) return NumericMatrix(0, 6);

    // Step 2-3: walk triplets, apply criteria, collect candidates.
    // Store as flat double vector (6 values per candidate).
    // Mirrors R: seq(1, length(zc)-3, by=2) -> 0-indexed: k = 0, 2, ..., nzc-4
    std::vector<double> flat;
    flat.reserve(nzc * 3);

    for (int k = 0; k < nzc - 3; k += 2) {
        int neg_start = zc[k];
        int neg_end   = zc[k + 1];
        int pos_end   = zc[k + 2];

        double dur_neg = (double)(neg_end - neg_start) / sr;
        double dur_pos = (double)(pos_end - neg_end)   / sr;

        if (dur_neg < dur_neg_min || dur_neg > dur_neg_max) continue;
        if (dur_pos < dur_pos_min || dur_pos > dur_pos_max) continue;

        double neg_peak = sig[neg_start];
        for (int j = neg_start + 1; j <= neg_end; j++)
            if (sig[j] < neg_peak) neg_peak = sig[j];

        double pos_peak = sig[neg_end];
        for (int j = neg_end + 1; j <= pos_end; j++)
            if (sig[j] > pos_peak) pos_peak = sig[j];

        double ptp = pos_peak - neg_peak;
        if (ptp < amp_min || ptp > amp_max) continue;

        // Store 1-indexed positions (matching R's which() convention)
        flat.push_back((double)(neg_start + 1));
        flat.push_back((double)(neg_end   + 1));
        flat.push_back((double)(pos_end   + 1));
        flat.push_back(neg_peak);
        flat.push_back(pos_peak);
        flat.push_back(ptp);
    }

    int nc = (int)flat.size() / 6;
    if (nc == 0) return NumericMatrix(0, 6);

    NumericMatrix out(nc, 6);
    for (int i = 0; i < nc; i++) {
        for (int j = 0; j < 6; j++) out(i, j) = flat[i * 6 + j];
    }

    Rcpp::colnames(out) = Rcpp::CharacterVector::create(
        "neg_start", "neg_end", "pos_end", "neg_peak", "pos_peak", "ptp"
    );
    return out;
}
