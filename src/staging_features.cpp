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
// [[Rcpp::export]]
NumericVector resample_poly_cpp(NumericVector x, NumericVector h, int up) {
    int N   = x.size();
    int nh  = h.size();
    int n_pre     = nh / 2;
    int max_taps  = (nh + up - 1) / up;
    int out_len   = N * up;

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
    for (int i = 0; i < out_len; i++) {
        int k         = i + n_pre;
        int phase     = k % up;
        int input_idx = k / up;
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
// [[Rcpp::export]]
double perm_entropy_cpp(NumericVector x, int order = 3, int delay = 1) {
    int N = x.size();
    int n = N - (order - 1) * delay;
    if (n <= 0) return NA_REAL;

    std::vector<int> fact(order, 1);
    for (int i = 1; i < order; i++) fact[i] = fact[i - 1] * i;

    int n_perms = fact[order - 1] * order;
    std::vector<int> counts(n_perms, 0);
    std::vector<int> argsort(order);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < order; j++) argsort[j] = j;
        std::stable_sort(argsort.begin(), argsort.end(),
            [&](int a, int b) {
                return x[i + a * delay] < x[i + b * delay];
            });
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

    double h = 0.0;
    for (int k = 0; k < n_perms; k++) {
        if (counts[k] > 0) {
            double p = (double)counts[k] / (double)n;
            h -= p * std::log(p);
        }
    }
    double log_n_perms = 0.0;
    for (int i = 1; i <= order; i++) log_n_perms += std::log((double)i);
    return h / log_n_perms;
}

// ── higuchi_fd_cpp ────────────────────────────────────────────────────────────
// Matches antropy.higuchi_fd(x, kmax=10) and R's .higuchi_fd.
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
            if (Lmk > 0.0) { Lk_sum += Lmk; valid_cnt++; }
        }
        if (valid_cnt > 0) {
            Lk[k - 1]       = Lk_sum / (double)valid_cnt;
            Lk_valid[k - 1] = true;
        }
    }

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
        sum_x  += lx[i]; sum_y  += ly[i];
        sum_xx += lx[i] * lx[i]; sum_xy += lx[i] * ly[i];
    }
    return (n_pts * sum_xy - sum_x * sum_y) /
           (n_pts * sum_xx - sum_x * sum_x);
}

// ── roll_triang_mean_cpp ──────────────────────────────────────────────────────
// Matches pandas rolling(window=k, center=True, min_periods=1, win_type='triang').mean()
//
// [[Rcpp::export]]
NumericVector roll_triang_mean_cpp(NumericVector x, int k = 15) {
    int n    = x.size();
    int half = (k - 1) / 2;
    std::vector<double> w_full(k);
    for (int j = 0; j <= half; j++)     w_full[j] = (double)(j + 1);
    for (int j = half + 1; j < k; j++) w_full[j] = (double)(k - j);

    NumericVector result(n);
    for (int i = 0; i < n; i++) {
        int i_start = std::max(0, i - half);
        int i_end   = std::min(n - 1, i + half);
        int w_start = i_start - (i - half);
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
// nzc_p = sign changes in consecutive first differences (local extrema count).
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
// hmob = sqrt(var(diff(x)) / var(x)), hcomp = sqrt(var(diff2(x)) / var(diff(x))) / hmob
// Single O(N) pass; Bessel-corrected variance to match R's var().
//
// [[Rcpp::export]]
NumericVector hjorth_cpp(NumericVector x) {
    int N = x.size();
    NumericVector na_out = NumericVector::create(
        Named("hmob") = NA_REAL, Named("hcomp") = NA_REAL
    );
    if (N < 3) return na_out;

    double sum_x  = 0, sum_x2  = 0;
    double sum_d1 = 0, sum_d12 = 0;
    double sum_d2 = 0, sum_d22 = 0;
    for (int i = 0; i < N; i++) {
        sum_x  += x[i];
        sum_x2 += x[i] * x[i];
        if (i < N - 1) {
            double d1 = x[i + 1] - x[i];
            sum_d1  += d1; sum_d12 += d1 * d1;
        }
        if (i < N - 2) {
            double d2 = x[i + 2] - 2.0 * x[i + 1] + x[i];
            sum_d2  += d2; sum_d22 += d2 * d2;
        }
    }
    int Nd1 = N - 1, Nd2 = N - 2;
    double var_x  = (sum_x2  - sum_x  * sum_x  / (double)N)   / (double)(N   - 1);
    double var_d1 = (sum_d12 - sum_d1 * sum_d1 / (double)Nd1) / (double)(Nd1 - 1);
    double var_d2 = (sum_d22 - sum_d2 * sum_d2 / (double)Nd2) / (double)(Nd2 - 1);

    double eps   = std::numeric_limits<double>::epsilon();
    double hmob  = std::sqrt(var_d1 / (var_x  + eps));
    double hcomp = std::sqrt(var_d2 / (var_d1 + eps)) / (hmob + eps);
    return NumericVector::create(Named("hmob") = hmob, Named("hcomp") = hcomp);
}

// ── stat_features_cpp ─────────────────────────────────────────────────────────
// std (Bessel), IQR (type-7), skewness mean(z^3), excess kurtosis mean(z^4)-3.
//
// [[Rcpp::export]]
NumericVector stat_features_cpp(NumericVector x) {
    int N = x.size();
    double eps = std::numeric_limits<double>::epsilon();

    double sum = 0.0;
    for (int i = 0; i < N; i++) sum += x[i];
    double mu = sum / (double)N;

    double sum2 = 0.0;
    for (int i = 0; i < N; i++) { double d = x[i] - mu; sum2 += d * d; }
    double var = (N > 1) ? sum2 / (double)(N - 1) : 0.0;
    double sig = std::sqrt(var);

    if (sig < eps) {
        return NumericVector::create(
            Named("std") = 0.0, Named("iqr") = 0.0,
            Named("skew") = 0.0, Named("kurt") = 0.0
        );
    }

    double sum3 = 0.0, sum4 = 0.0;
    for (int i = 0; i < N; i++) {
        double z = (x[i] - mu) / sig, z2 = z * z;
        sum3 += z2 * z; sum4 += z2 * z2;
    }

    std::vector<double> xs(x.begin(), x.end());
    std::sort(xs.begin(), xs.end());
    auto quant7 = [&](double p) -> double {
        double h = (double)(N - 1) * p;
        int lo = (int)std::floor(h), hi = lo + 1;
        double frac = h - (double)lo;
        if (hi >= N) return xs[N - 1];
        return xs[lo] + frac * (xs[hi] - xs[lo]);
    };

    return NumericVector::create(
        Named("std")  = sig,
        Named("iqr")  = quant7(0.75) - quant7(0.25),
        Named("skew") = sum3 / (double)N,
        Named("kurt") = sum4 / (double)N - 3.0
    );
}

// ── rowmedian_cpp ─────────────────────────────────────────────────────────────
// Row-wise median of a numeric matrix.
// Replaces apply(pgrams, 1L, stats::median) in .welch_median_psd —
// avoids 251 R dispatch calls per epoch.
//
// [[Rcpp::export]]
NumericVector rowmedian_cpp(NumericMatrix x) {
    int nrow = x.nrow();
    int ncol = x.ncol();
    NumericVector out(nrow);
    std::vector<double> buf(ncol);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) buf[j] = x(i, j);
        std::sort(buf.begin(), buf.end());
        out[i] = (ncol % 2 == 1) ? buf[ncol / 2]
                                  : (buf[ncol / 2 - 1] + buf[ncol / 2]) / 2.0;
    }
    return out;
}

// ── roll_right_mean_cpp ───────────────────────────────────────────────────────
// Right-aligned rolling mean with partial windows at the start.
// Matches zoo::rollapply(x, k, mean, fill=NA, partial=TRUE, align="right").
//
// [[Rcpp::export]]
NumericVector roll_right_mean_cpp(NumericVector x, int k = 4) {
    int n = x.size();
    NumericVector out(n);
    for (int i = 0; i < n; i++) {
        int start = std::max(0, i - k + 1);
        int count = i - start + 1;
        double sum = 0.0;
        for (int j = start; j <= i; j++) sum += x[j];
        out[i] = sum / (double)count;
    }
    return out;
}

// ── robust_scale_cpp ──────────────────────────────────────────────────────────
// (x - median) / (quantile(q_high) - quantile(q_low) + 1e-10).
// NA-aware: non-NA values used for stats; NA positions preserved in output.
// Type-7 quantile matching R's default quantile().
//
// [[Rcpp::export]]
NumericVector robust_scale_cpp(NumericVector x,
                                double q_low  = 0.05,
                                double q_high = 0.95) {
    int n_total = x.size();
    std::vector<double> xs;
    xs.reserve(n_total);
    for (int i = 0; i < n_total; i++) {
        if (!ISNAN(x[i])) xs.push_back(x[i]);
    }
    int n = xs.size();
    if (n == 0) return NumericVector(n_total, NA_REAL);

    std::sort(xs.begin(), xs.end());
    double med = (n % 2 == 1) ? xs[n / 2]
                               : (xs[n / 2 - 1] + xs[n / 2]) / 2.0;

    auto quant7 = [&](double p) -> double {
        double h = (double)(n - 1) * p;
        int lo = (int)std::floor(h), hi = lo + 1;
        double frac = h - (double)lo;
        if (hi >= n) return xs[n - 1];
        return xs[lo] + frac * (xs[hi] - xs[lo]);
    };
    double scale = quant7(q_high) - quant7(q_low) + 1e-10;

    NumericVector out(n_total);
    for (int i = 0; i < n_total; i++) {
        out[i] = ISNAN(x[i]) ? NA_REAL : (x[i] - med) / scale;
    }
    return out;
}
