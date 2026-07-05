// src/staging_features.cpp
//
// C++ implementations of the per-epoch feature functions that are bottlenecks
// in pure R: perm_entropy, higuchi_fd, roll_triang_mean, and resample_poly.
//
// Each function mirrors its R/Python counterpart exactly — same algorithm, same
// numerical output — but at C++ speed.  The R wrappers in staging_features.R
// call these directly; the pure-R fallbacks are kept for reference.
//
// Build: Rcpp::compileAttributes()  then  devtools::load_all()

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// ── resample_poly_cpp ────────────────────────────────────────────────────────────────────
// Polyphase integer upsampler matching scipy.signal.resample_poly(x, up, 1).
//
// scipy.signal.resample_poly calls upfirdn(h, x, up, 1) then trims the output:
//   full_out = upfirdn(h, x, up, 1)       # length = N*up + nh - 1
//   out      = full_out[n_pre : n_pre + N*up]  where n_pre = nh // 2
//
// The polyphase decomposition of h into `up` branches:
//   poly[p][j] = h[p + j*up]   for p = 0..up-1, j = 0..ceil(nh/up)-1
//
// upfirdn output at position k:
//   y[k] = sum_j h[k%up + j*up] * x[k//up - j]   (zero-padding for out-of-bounds x)
//
// By computing y[n_pre .. n_pre+N*up-1] directly, we exactly replicate the
// scipy/upfirdn output without materialising the zero-inserted signal.
//
// [[Rcpp::export]]
NumericVector resample_poly_cpp(NumericVector x, NumericVector h, int up) {
    int N   = x.size();
    int nh  = h.size();            // 2001
    int n_pre     = nh / 2;        // = 1000  (filter half-delay in output samples)
    int max_taps  = (nh + up - 1) / up;   // = ceil(2001/100) = 21
    int out_len   = N * up;

    // Build polyphase matrix: poly[p][j] = h[p + j*up]
    std::vector<std::vector<double>> poly(
        up, std::vector<double>(max_taps, 0.0)
    );
    for (int p = 0; p < up; p++) {
        for (int j = 0; j < max_taps; j++) {
            int idx = p + j * up;
            if (idx < nh) poly[p][j] = h[idx];
        }
    }

    NumericVector out(out_len);

    // Compute y[k] for k = n_pre .. n_pre + out_len - 1
    // (equivalent to scipy's trim: full[n_pre : n_pre + N*up])
    for (int i = 0; i < out_len; i++) {
        int k         = i + n_pre;   // position in full upfirdn output
        int phase     = k % up;
        int input_idx = k / up;      // floor(k / up)

        double sum = 0.0;
        for (int j = 0; j < max_taps; j++) {
            int in_j = input_idx - j;
            if (in_j >= 0 && in_j < N) sum += poly[phase][j] * x[in_j];
        }
        out[i] = sum;
    }
    return out;
}

// ── perm_entropy_cpp ─────────────────────────────────────────────────────────
// Matches antropy.perm_entropy(x, order=3, delay=1, normalize=True).
//
// Instead of R's string-based pattern hashing (apply / paste / table),
// we encode each permutation as an integer via the Lehmer (factoriadic) code:
//   for position i in the argsort, count how many elements to its right
//   have a smaller value.  This gives a unique index in [0, order! - 1].
//
// [[Rcpp::export]]
double perm_entropy_cpp(NumericVector x, int order = 3, int delay = 1) {
    int N = x.size();
    int n = N - (order - 1) * delay;
    if (n <= 0) return NA_REAL;

    // Precompute factorials for Lehmer encoding
    std::vector<int> fact(order, 1);
    for (int i = 1; i < order; i++) fact[i] = fact[i - 1] * i;

    int n_perms = fact[order - 1] * order;   // order!
    std::vector<int> counts(n_perms, 0);

    std::vector<int> argsort(order);

    for (int i = 0; i < n; i++) {
        // argsort for window [x[i], x[i+delay], ..., x[i+(order-1)*delay]]
        for (int j = 0; j < order; j++) argsort[j] = j;

        // Stable sort ascending — matches R's order() which is also stable
        std::stable_sort(argsort.begin(), argsort.end(),
            [&](int a, int b) {
                return x[i + a * delay] < x[i + b * delay];
            });

        // Lehmer code: for each position j, count elements to the right
        // with smaller value in argsort (= factoriadic encoding)
        int hash = 0;
        for (int j = 0; j < order - 1; j++) {
            int cnt = 0;
            for (int k = j + 1; k < order; k++) {
                if (argsort[k] < argsort[j]) cnt++;
            }
            hash += cnt * fact[order - 1 - j];
        }
        counts[hash]++;
    }

    // Shannon entropy over observed pattern frequencies
    double h = 0.0;
    for (int k = 0; k < n_perms; k++) {
        if (counts[k] > 0) {
            double p = (double)counts[k] / (double)n;
            h -= p * std::log(p);
        }
    }

    // Normalise by log(order!) — matches antropy normalize=True and R version
    double log_n_perms = 0.0;
    for (int i = 1; i <= order; i++) log_n_perms += std::log((double)i);

    return h / log_n_perms;
}

// ── higuchi_fd_cpp ────────────────────────────────────────────────────────────
// Matches antropy.higuchi_fd(x, kmax=10) and R's .higuchi_fd.
//
// Direct C++ translation of the nested loop; the OLS slope is computed
// via the standard closed-form formula to avoid any R lm() overhead.
//
// [[Rcpp::export]]
double higuchi_fd_cpp(NumericVector x, int kmax = 10) {
    int N = x.size();

    std::vector<double> Lk(kmax, 0.0);
    std::vector<bool>   Lk_valid(kmax, false);

    for (int k = 1; k <= kmax; k++) {
        double Lk_sum    = 0.0;
        int    valid_cnt = 0;

        for (int m = 1; m <= k; m++) {
            // R: idx <- seq.int(m, N, by = k)  (1-based)
            // 0-based: m-1, m-1+k, m-1+2k, ...
            // km = floor((N - m) / k) = number of step differences
            int km = (N - m) / k;
            if (km < 1) continue;

            double sum_abs = 0.0;
            for (int t = 0; t < km; t++) {
                int idx_next = (m - 1) + (t + 1) * k;
                int idx_curr = (m - 1) +  t      * k;
                sum_abs += std::abs(x[idx_next] - x[idx_curr]);
            }

            double Lmk = sum_abs * (double)(N - 1) /
                         ((double)km * (double)k * (double)k);

            if (Lmk > 0.0) {
                Lk_sum += Lmk;
                valid_cnt++;
            }
        }

        if (valid_cnt > 0) {
            Lk[k - 1]      = Lk_sum / (double)valid_cnt;
            Lk_valid[k - 1] = true;
        }
    }

    // Collect valid (k, Lk) pairs for OLS on log-log scale
    std::vector<double> lx, ly;
    for (int k = 1; k <= kmax; k++) {
        if (Lk_valid[k - 1] && Lk[k - 1] > 0.0) {
            lx.push_back(std::log(1.0 / (double)k));
            ly.push_back(std::log(Lk[k - 1]));
        }
    }

    if ((int)lx.size() < 2) return NA_REAL;

    double n_pts  = (double)lx.size();
    double sum_x  = 0, sum_y  = 0, sum_xx = 0, sum_xy = 0;
    for (int i = 0; i < (int)lx.size(); i++) {
        sum_x  += lx[i];
        sum_y  += ly[i];
        sum_xx += lx[i] * lx[i];
        sum_xy += lx[i] * ly[i];
    }

    return (n_pts * sum_xy - sum_x * sum_y) /
           (n_pts * sum_xx - sum_x * sum_x);
}

// ── roll_triang_mean_cpp ──────────────────────────────────────────────────────
// Matches pandas rolling(window=k, center=True, min_periods=1, win_type='triang').mean()
// and R's .roll_triang_mean.
//
// Triangle weights for a full window of k:
//   w[j] = j+1        for j in 0 .. half
//   w[j] = k - j      for j in half+1 .. k-1
// Edge windows use the corresponding slice of w, then normalise.
//
// [[Rcpp::export]]
NumericVector roll_triang_mean_cpp(NumericVector x, int k = 15) {
    int n    = x.size();
    int half = (k - 1) / 2;

    // Full triangular weight vector (length k)
    std::vector<double> w_full(k);
    for (int j = 0; j <= half; j++)     w_full[j] = (double)(j + 1);
    for (int j = half + 1; j < k; j++) w_full[j] = (double)(k - j);

    NumericVector result(n);

    for (int i = 0; i < n; i++) {
        int i_start = std::max(0, i - half);
        int i_end   = std::min(n - 1, i + half);
        int w_start = i_start - (i - half);   // offset into w_full

        double wsum = 0.0, w_total = 0.0;
        for (int j = i_start; j <= i_end; j++) {
            double w  = w_full[w_start + (j - i_start)];
            wsum    += x[j] * w;
            w_total += w;
        }
        result[i] = wsum / w_total;
    }
    return result;
}
