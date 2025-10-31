# install.packages("multitaper")  # uncomment if needed
library(multitaper)
library(tidyverse)

options(scipen = 6) # print fewer numbers in scientific notation (more plain decimals)
set.seed(123) # reproducibility for the toy simulation


data.names <-
  list.files("../Data", full.names = T, pattern = "*.csv$")
data.files <- read_csv(data.names)

data.aggr <-
  data.files |>
  filter(Trial_type == "experiment", Rt > 0.0) |> # Experiment and no timeout
  dplyr::select(PART_ID, CSI, Corr) |>
  mutate(t = CSI * 16.6) |>
  group_by(PART_ID, t) |>
  summarise(mean_corr = 100.0 * mean(Corr))


# =====================================
# 1) Toy data: MANY short time series
# =====================================
# Goal: simulate M short time series, some with a narrowband sinusoid + noise, others pure noise.
# We'll later test each series for a sinusoidal (line) component using the Thomson multitaper harmonic F-test.
# References: Thomson (1982) Proceedings of the IEEE; Percival & Walden (1993) CUP.

M <- 30 # number of independent time series (e.g., subjects, trials, sensors). Use 150 in your real case.
N <- 256 # samples per series (shorter segments benefit from multitaper for variance reduction)
fs <- 1 # sampling rate in Hz (so Nyquist is 0.5 Hz)
dt <- 1 / fs
t <- seq(0, (N - 1) / fs, by = dt) # time vector

# Mix: ~60% with a random sinusoid + noise, ~40% pure noise
# This creates ground truth to evaluate detection performance (TP, FP, etc.) in this toy example.
has_signal <- rbinom(M, 1, 0.6) == 1
series <- vector("list", M)
true_freqs <- rep(NA_real_, M) # store the true sinusoid frequency for signal series (NA otherwise)

for (i in 1:M) {
  if (has_signal[i]) {
    # draw a random sinusoid frequency away from 0 (DC) and Nyquist to avoid edge issues
    f0 <- runif(1, 0.08, 0.45)
    A <- runif(1, 0.8, 1.6) # sinusoid amplitude
    phi <- runif(1, 0, 2 * pi) # random phase
    # signal + white Gaussian noise (sd = 1.0)
    x <- A * sin(2 * pi * f0 * t + phi) + rnorm(N, 0, 1.0)
    true_freqs[i] <- f0
  } else {
    # pure noise series (null hypothesis true)
    x <- rnorm(N, 0, 1.0)
  }
  series[[i]] <- x
}

# ==========================================================
# 2) MTM parameters and helper (Thomson 1982 harmonic F-test)
# ==========================================================
# Multitaper (MTM) core ideas:
# - Use K orthogonal DPSS tapers to get K tapered spectra; average to reduce variance without excessive smoothing.
# - Time–bandwidth product NW controls spectral leakage and resolution. Half-bandwidth W = NW / (N*dt).
# - For detecting a sinusoidal "line" (discrete spectral component), Thomson (1982) derived an F-test whose null is
#   "no line at that frequency" (i.e., locally broadband) and alternative is "line present".
#
# Parameter choices:
# - NW = 3 is a common choice for short records; K = 2*NW - 1 tapers is a widely used rule of thumb.
#   Tradeoff: larger NW increases smoothing (reduces variance) but broadens spectral features.
#   Here, W = 3 / (256 * 1) ≈ 0.0117 Hz half-bandwidth.
# References:
# - Thomson, D. J. (1982). Spectrum estimation and harmonic analysis. Proc. IEEE. https://doi.org/10.1109/PROC.1982.12433
# - Percival, D. B., & Walden, A. T. (1993). Spectral Analysis for Physical Applications. Cambridge University Press.

NW <- 3
K <- 2 * NW - 1 # number of tapers used (here 5), commonly used setting
alpha <- 0.05 # significance level for reporting (per-series and global summary)

# Helper to robustly extract the F-test vector across multitaper versions/structures
get_F <- function(obj) {
  # The spec.mtm object can store F-test results in different nested locations depending on package/version.
  if (!is.null(obj$mtm$Ftest)) {
    return(obj$mtm$Ftest)
  }
  if (!is.null(obj$Ftest)) {
    return(obj$Ftest)
  }
  if (!is.null(obj$mtm$mtm$Ftest)) {
    return(obj$mtm$mtm$Ftest)
  }
  stop("Could not locate F-test values in the spec.mtm object. Inspect with str(obj).")
}

# Degrees of freedom for Thomson harmonic F-test:
# Under H0 (no line), at each frequency the F-statistic approximately follows F_{2, 2K-2}
# (2 numerator df from the 2 parameters of the sinusoid; denominator df from K tapers)
df1 <- 2
df2 <- 2 * K - 2 # with K=5, df2=8

# ===================================
# 3) Compute MTM + F-test per series
# ===================================
# For each series:
# - Compute multitaper spectrum with Thomson's harmonic F-test enabled.
# - Convert F to per-frequency p-values using the F distribution.
# Notes:
# - Frequencies within ~W are not independent due to taper bandwidth; we handle this later via a conservative Bonferroni.
# - spec.mtm frequency grid is in cycles per unit time (Hz here because dt is in seconds).
mtm_list <- vector("list", M)
freq <- NULL

for (i in 1:M) {
  mtm <- spec.mtm(ts(series[[i]], deltat = dt),
    nw = NW, k = K, Ftest = TRUE, plot = FALSE,
    deltat = dt, dtUnits = "s"
  )
  Fvals <- get_F(mtm)
  pvals <- 1 - pf(Fvals, df1, df2) # per-frequency p-values for H0: no line at that frequency
  mtm_list[[i]] <- list(freq = mtm$freq, F = Fvals, p = pvals)
  if (is.null(freq)) freq <- mtm$freq # store common frequency grid from the first series
}

# Optional: restrict frequency band (avoid DC and very high end)
# Rationale: DC (0 Hz) often includes mean/trend; near Nyquist can have edge artifacts.
keep <- which(freq > 0 & freq < (fs / 2 - 1e-6))
m <- length(keep) # number of tested frequencies per series (used for Bonferroni within-series)

# =====================================================================
# 4) Reduce each series to one p-value (Bonferroni across frequencies)
# =====================================================================
# We scan m frequencies in each series and want a single per-series decision about "any line present".
# Strategy:
# - Take the minimum per-frequency p-value and apply a Bonferroni correction: p_bonf = min(1, m * min_p_raw).
# - This controls family-wise error rate (FWER) within each series under independence. With dependence (here present),
#   Bonferroni is conservative (erring on the side of fewer false positives).
# Alternatives (not implemented here):
# - Use the maximum F across frequencies and simulate its null via parametric/bootstrap to respect dependence.
# - Use Simes or Holm across frequencies within series.
p_series <- numeric(M) # Bonferroni-corrected per-series p-value
min_p_raw <- numeric(M) # minimal unadjusted per-frequency p-value
peak_freq <- numeric(M) # frequency at which min p (and max F) occurs
peak_F <- numeric(M) # corresponding F-statistic
peak_index <- integer(M) # index of that frequency

for (i in 1:M) {
  p_i <- mtm_list[[i]]$p[keep]
  F_i <- mtm_list[[i]]$F[keep]
  idx_min <- which.min(p_i) # frequency with the strongest evidence for a line
  pmin <- p_i[idx_min]
  # Bonferroni across m tested frequencies within this series
  p_bonf <- min(1, m * pmin)

  p_series[i] <- p_bonf
  min_p_raw[i] <- pmin
  peak_index[i] <- keep[idx_min]
  peak_freq[i] <- freq[keep][idx_min]
  peak_F[i] <- F_i[idx_min]
}

# ===================================================
# 5) Combine per-series p-values into ONE decision
#     Fisher's method (global p-value)
# ===================================================
# We now ask: "Is there evidence that ANY of the M series contain a line?"
# Combine the M per-series (FWER-adjusted) p-values via Fisher’s method:
#   T = -2 * sum(log p_i) ~ Chi^2_{2M} under independence.
# Caveat: The independence assumption across series is often reasonable (e.g., different trials/subjects),
# but if series are dependent, Fisher can be anti-conservative. There are robust alternatives (e.g., Cauchy combination).
# References:
# - Fisher, R.A. (1932). Statistical Methods for Research Workers.
# - Brown, M. B. (1975). A Method for Combining Non-Independent, One-Sided Tests. Biometrics. https://doi.org/10.2307/2529716
# - Liu, Y., & Xie, J. (2020). Cauchy Combination Test. JASA. https://doi.org/10.1080/01621459.2020.1758115
T_fisher <- -2 * sum(log(p_series))
p_global <- 1 - pchisq(T_fisher, df = 2 * M)

# Optional robust alternative (Cauchy combination) that tolerates arbitrary dependency among p-values:
# T_c  <- mean(tan((0.5 - p_series)*pi))
# p_c  <- 0.5 - atan(T_c)/pi
# (Not used in the printed report, but provided for reference.)

# ===============================================
# 6) Formatted report printed to the console
# ===============================================
# Pretty-print a compact report with settings, global decision, and per-series summary.
line <- function(char = "=", n = 70) cat(paste0(paste(rep(char, n), collapse = ""), "\n"))
fmt <- function(x, d = 3) formatC(x, digits = d, format = "fg", flag = "#")

line("=")
cat("Multitaper Harmonic F-test: Group Decision Report\n")
line("=")
cat("Data settings:\n")
cat(sprintf("  Series count (M): %d\n", M))
cat(sprintf("  Samples per series (N): %d\n", N))
cat(sprintf("  Sampling rate (fs): %s Hz (dt = %s s)\n", fmt(fs, 6), fmt(dt, 6)))
cat(sprintf("  Frequency bins tested per series (m): %d\n", m))
cat(sprintf(
  "  MTM parameters: NW = %s, K = %d (df1=%d, df2=%d)\n",
  fmt(NW, 3), K, df1, df2
))
line("-")
cat("Global test (Fisher over per-series Bonferroni p-values):\n")
cat(sprintf(
  "  Test statistic T = %s, df = %d, p_global = %s\n",
  fmt(T_fisher, 4), 2 * M, fmt(p_global, 4)
))
cat(sprintf(
  "  Decision at alpha = %.2f: %s\n",
  alpha,
  ifelse(p_global < alpha,
    "REJECT H0 (oscillations present in the collection)",
    "FAIL TO REJECT H0 (no evidence for oscillations)"
  )
))
line("-")

# Per-series summary table
# 'sig_at_alpha' reflects whether the series-level Bonferroni p-value <= alpha.
# Note: This is a within-series multiple-comparison correction. It does not adjust across M series (that is the job of the global test).
sig_series_flag <- p_series <= alpha
summary_df <- data.frame(
  series_id    = 1:M,
  has_signal   = ifelse(is.na(true_freqs), FALSE, TRUE), # ground truth available only in this toy setup
  peak_freq    = round(peak_freq, 4),
  peak_F       = round(peak_F, 3),
  min_p_raw    = signif(min_p_raw, 3),
  p_bonf       = signif(p_series, 3),
  sig_at_alpha = sig_series_flag
)

cat("Per-series summary (Bonferroni-corrected p over frequency search):\n")
print(summary_df, row.names = FALSE)
line("-")

# Counts and quick diagnostics (available in this toy setup)
# Compute confusion-matrix-style diagnostics leveraging ground truth:
# TP: detected and truly has a line; FP: detected but actually noise; TN, FN similarly.
if (any(!is.na(true_freqs))) {
  TP <- sum(sig_series_flag & has_signal)
  FP <- sum(sig_series_flag & !has_signal)
  TN <- sum(!sig_series_flag & !has_signal)
  FN <- sum(!sig_series_flag & has_signal)
  cat("Toy-simulation diagnostics (since ground truth is known here):\n")
  cat(sprintf(
    "  Significant series: %d/%d (TP=%d, FP=%d, TN=%d, FN=%d)\n",
    sum(sig_series_flag), M, TP, FP, TN, FN
  ))
  # Optional: frequency error for detected signals
  detected <- which(sig_series_flag & has_signal)
  if (length(detected) > 0) {
    freq_err <- abs(true_freqs[detected] - peak_freq[detected])
    cat(sprintf("  Median |freq_error| over detected signals: %s Hz\n", fmt(median(freq_err), 5)))
  }
  line("-")
}

# Top-10 strongest series by per-series evidence (smallest Bonferroni p)
# Useful to inspect which series drive the global rejection.
ord <- order(p_series)
kshow <- min(10, M)
cat(sprintf("Top %d series by strength (smallest per-series Bonferroni p):\n", kshow))
print(summary_df[ord[1:kshow], ], row.names = FALSE)
line("=")

# =================================
# 7) Visualizations (quick panels)
# =================================
# Four simple plots to sanity-check per-series results and the global decision:
# (a) Histogram of per-series Bonferroni p-values with alpha line.
# (b) QQ-plot vs Uniform(0,1): under global null, p-values ~ U(0,1); deviations indicate signal presence.
# (c) Scatter of peak frequency vs peak F; color red if series is significant at alpha.
# (d) For signal series only (toy), compare true vs detected peak frequency.

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

# (a) Histogram of per-series Bonferroni p-values
hist(p_series,
  breaks = 12, col = "gray85", border = "white",
  main = "Per-series Bonferroni-adjusted p-values", xlab = "p_i*"
)
abline(v = alpha, col = "red", lty = 2)
mtext(sprintf("Global p (Fisher) = %s", fmt(p_global, 4)), side = 3, line = 0.2, cex = 0.9)

# (b) QQ-plot vs Uniform(0,1)
# If most points fall below the diagonal, there is an excess of small p-values → evidence for signals.
pp <- sort(p_series)
uu <- (1:M) / (M + 1)
plot(uu, pp,
  pch = 19, cex = 0.8,
  xlab = "Theoretical U(0,1) quantiles", ylab = "Observed per-series p-values",
  main = "QQ plot of per-series p-values"
)
abline(0, 1, col = "gray50")

# (c) Peak frequency vs peak F (red if significant per-series)
# Stronger peaks (larger F) and mid-band frequencies often indicate clearer sinusoidal lines.
plot(peak_freq, peak_F,
  pch = 19, col = ifelse(sig_series_flag, "red", "black"),
  xlab = "Peak frequency (Hz)", ylab = "Peak F-statistic",
  main = "Strongest line per series"
)
legend("topright",
  legend = c("significant (per-series)", "non-significant"),
  pch = 19, col = c("red", "black"), bty = "n", cex = 0.9
)

# (d) True vs detected peak frequency (only for signal series in toy data)
# Ideally, points lie on the diagonal (accurate frequency estimation). Misses (non-significant) plotted in gray.
idx_sig <- which(!is.na(true_freqs))
if (length(idx_sig) > 0) {
  plot(true_freqs[idx_sig], peak_freq[idx_sig],
    xlab = "True frequency (Hz)", ylab = "Detected peak frequency (Hz)",
    main = "True vs detected peak (signal series only)",
    pch = 19, col = ifelse(sig_series_flag[idx_sig], "red", "gray50")
  )
  abline(0, 1, col = "gray60", lty = 2)
} else {
  plot.new()
  title("No ground truth in this dataset")
}

# Notes and practical tips:
# - Choice of NW/K: For N=256, NW in [2.5, 4] and K=2*NW-1 is common. If signals are very narrow, smaller NW improves resolution
#   but increases variance. Always check robustness with nearby settings.
# - Frequency search multiplicity: Bonferroni is conservative because frequency bins are dependent (~W correlation). If power is
#   critical, consider dependence-aware max-statistic calibration or within-series FDR methods.
# - Across-series multiplicity: Fisher’s method tests the global null (all series null) vs. "at least one non-null". If you instead
#   want to control FWER/FDR over the list of series called significant, apply Holm/BH to the per-series p-values.
# - Interpretability: The global p-value summarizes evidence that there exists oscillatory structure somewhere in the collection,
#   while per-series summaries tell you where and at what frequency peaks occur.
# - References:
#   * Thomson (1982): https://doi.org/10.1109/PROC.1982.12433
#   * Percival & Walden (1993): https://doi.org/10.1017/CBO9780511622762
#   * Fisher’s method: https://doi.org/10.2307/2342435 (historical); Brown (1975): https://doi.org/10.2307/2529716
#   * Cauchy combination test: https://doi.org/10.1080/01621459.2020.1758115
