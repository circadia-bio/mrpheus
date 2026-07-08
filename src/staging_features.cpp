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

// ── nzc_cpp ───────────────────────────────────────────────────────────────────
// Matches R's .nzc: sum(diff(sign(x)) != 0L)
// Counts the number of times the sign of x changes between consecutive samples
// (zero crossings, including crossings through zero).
//
// sign(v) is -1 / 0 / +1 using the C ternary trick: (v > 0) - (v < 0).
//
// [[Rcpp::export]]
int nzc_cpp(NumericVector x) {
    int N     = x.size();
    int count = 0;
    for (int i = 0; i < N - 1; i++) {
        int s1 = (x[i]     > 0.0) - (x[i]     < 0.0);
        int s2 = (x[i + 1] > 0.0) - (x[i + 1] < 0.0);
        if (s1 != s2) count++;
    }
    return count;
}

// ── petrosian_fd_cpp ──────────────────────────────────────────────────────────
// Matches antropy.petrosian_fd and R's .petrosian_fd.
//
// nzc_p = sign changes in *consecutive first differences* of x
//       = local extrema count (NOT zero crossings of x itself).
// R equivalent: dx <- diff(x); sum((dx[-length(dx)] * dx[-1]) < 0)
//
// [[Rcpp::export]]
double petrosian_fd_cpp(NumericVector x) {
    int N = x.size();
    if (N < 3) return NA_REAL;

    int nzc_p = 0;
    for (int i = 0; i < N - 2; i++) {
        double d_curr = x[i + 1] - x[i];
        double d_next = x[i + 2] - x[i + 1];
        if (d_curr * d_next < 0.0) nzc_p++;
    }

    double logN = std::log10((double)N);
    return logN / (logN + std::log10((double)N /
                                     ((double)N + 0.4 * (double)nzc_p)));
}

// ── hjorth_cpp ────────────────────────────────────────────────────────────────
// Matches R's .hjorth:
//   hmob  = sqrt(var(diff(x))        / var(x))
//   hcomp = sqrt(var(diff(diff(x))) / var(diff(x))) / hmob
//
// Bessel-corrected variance (divides by N-1) to match R's var().
// Single O(N) pass avoids materialising diff() and diff(diff()) vectors.
// Returns a named NumericVector c(hmob = ..., hcomp = ...).
//
// [[Rcpp::export]]
NumericVector hjorth_cpp(NumericVector x) {
    int N = x.size();
    NumericVector na_out = NumericVector::create(
        Named("hmob") = NA_REAL, Named("hcomp") = NA_REAL
    );
    if (N < 3) return na_out;

    // Accumulate sums for x, d1=diff(x), d2=diff(diff(x)) in one pass
    double sum_x  = 0, sum_x2  = 0;
    double sum_d1 = 0, sum_d12 = 0;
    double sum_d2 = 0, sum_d22 = 0;

    for (int i = 0; i < N; i++) {
        sum_x  += x[i];
        sum_x2 += x[i] * x[i];
        if (i < N - 1) {
            double d1 = x[i + 1] - x[i];
            sum_d1  += d1;
            sum_d12 += d1 * d1;
        }
        if (i < N - 2) {
            double d2 = x[i + 2] - 2.0 * x[i + 1] + x[i];
            sum_d2  += d2;
            sum_d22 += d2 * d2;
        }
    }

    int Nd1 = N - 1, Nd2 = N - 2;

    // Bessel-corrected variance: (sum_sq - sum^2/n) / (n-1)
    double var_x  = (sum_x2  - sum_x  * sum_x  / (double)N)   / (double)(N   - 1);
    double var_d1 = (sum_d12 - sum_d1 * sum_d1 / (double)Nd1) / (double)(Nd1 - 1);
    double var_d2 = (sum_d22 - sum_d2 * sum_d2 / (double)Nd2) / (double)(Nd2 - 1);

    double eps   = std::numeric_limits<double>::epsilon();
    double hmob  = std::sqrt(var_d1 / (var_x  + eps));
    double hcomp = std::sqrt(var_d2 / (var_d1 + eps)) / (hmob + eps);

    return NumericVector::create(Named("hmob") = hmob, Named("hcomp") = hcomp);
}

// ── stat_features_cpp ─────────────────────────────────────────────────────────
// Matches R's .stat_features: std, IQR (type-7 quantile), skewness, kurtosis.
// Returns a named NumericVector c(std = ..., iqr = ..., skew = ..., kurt = ...).
//
// std  : Bessel-corrected (matches R's sd()).
// IQR  : R type-7 linear interpolation on sorted copy — h = (n-1)*p (0-indexed).
// skew : mean(z^3) where z = (x - mu) / sd.
// kurt : mean(z^4) - 3  (excess kurtosis).
//
// [[Rcpp::export]]
NumericVector stat_features_cpp(NumericVector x) {
    int N = x.size();
    double eps = std::numeric_limits<double>::epsilon();

    // Mean (single pass)
    double sum = 0.0;
    for (int i = 0; i < N; i++) sum += x[i];
    double mu = sum / (double)N;

    // Bessel-corrected variance (two-pass for numerical stability)
    double sum2 = 0.0;
    for (int i = 0; i < N; i++) {
        double d = x[i] - mu;
        sum2 += d * d;
    }
    double var = (N > 1) ? sum2 / (double)(N - 1) : 0.0;
    double sig = std::sqrt(var);

    if (sig < eps) {
        return NumericVector::create(
            Named("std")  = 0.0, Named("iqr")  = 0.0,
            Named("skew") = 0.0, Named("kurt") = 0.0
        );
    }

    // Skewness and excess kurtosis (normalised moments)
    double sum3 = 0.0, sum4 = 0.0;
    for (int i = 0; i < N; i++) {
        double z  = (x[i] - mu) / sig;
        double z2 = z * z;
        sum3 += z2 * z;
        sum4 += z2 * z2;
    }
    double skew = sum3 / (double)N;
    double kurt = sum4 / (double)N - 3.0;

    // IQR via R type-7 quantile on a sorted copy
    // h = (n-1)*p  (0-indexed float position), then linearly interpolate.
    std::vector<double> xs(x.begin(), x.end());
    std::sort(xs.begin(), xs.end());

    auto quant7 = [&](double p) -> double {
        double h    = (double)(N - 1) * p;
        int    lo   = (int)std::floor(h);
        int    hi   = lo + 1;
        double frac = h - (double)lo;
        if (hi >= N) return xs[N - 1];
        return xs[lo] + frac * (xs[hi] - xs[lo]);
    };

    double iqr = quant7(0.75) - quant7(0.25);

    return NumericVector::create(
        Named("std")  = sig,
        Named("iqr")  = iqr,
        Named("skew") = skew,
        Named("kurt") = kurt
    );
}
