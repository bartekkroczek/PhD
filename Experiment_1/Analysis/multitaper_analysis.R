library(dplyr)
library(tibble)
library(tidyverse)

# ============
# Parameters
# ============
FS <- 60
DT <- 1 / FS
NW <- 3
K <- 2 * NW - 1
alpha <- 0.05
fmin <- 0.0
fmax <- NULL # if NULL, uses fs/2 - 1e-6 per series

# F-test degrees of freedom (Thomson 1982)
df1 <- 2
df2 <- 2 * K - 2

# ============
# Helpers
# ============
line <- function(char = "=", n = 70) cat(paste0(paste(rep(char, n), collapse = ""), "\n"))
fmt <- function(x, d = 3) formatC(x, digits = d, format = "fg", flag = "#")

get_F <- function(obj) {
  if (!is.null(obj$mtm$Ftest)) {
    return(obj$mtm$Ftest)
  }
  if (!is.null(obj$Ftest)) {
    return(obj$Ftest)
  }
  if (!is.null(obj$mtm$mtm$Ftest)) {
    return(obj$mtm$mtm$Ftest)
  }
  stop("Could not locate F-test values in spec.mtm object. Inspect with str(obj).")
}

# ============
# Data prep
# ============

data.names <-
  list.files("../Data", full.names = T, pattern = "*.csv$")
data.files <- read_csv(data.names)

data.aggr <-
  data.files |>
  filter(Trial_type == "experiment", Rt > 0.0) |> # Experiment and no timeout
  dplyr::select(PART_ID, CSI, Corr) |>
  mutate(t = CSI * DT) |>
  group_by(PART_ID, t) |>
  summarise(mean_corr = 100.0 * mean(Corr))

stopifnot(all(c("PART_ID", "t", "mean_corr") %in% names(data.aggr)))

data_clean <- data.aggr %>%
  filter(is.finite(t), is.finite(mean_corr)) %>%
  arrange(PART_ID, t)

series_info <- data_clean %>%
  group_by(PART_ID) %>%
  summarize(N = n(), .groups = "drop") %>%
  mutate(dt = DT, fs = FS)

# Construct vectors per series
series_list <- data_clean %>%
  group_by(PART_ID) %>%
  arrange(t, .by_group = TRUE) %>%
  summarize(
    x = list(mean_corr),
    .groups = "drop"
  ) %>%
  left_join(series_info, by = "PART_ID")

M <- nrow(series_list)
if (M == 0L) stop("No series found after cleaning.")
if (any(!is.finite(series_list$dt) | series_list$dt <= 0)) {
  stop("Non-finite or non-positive dt inferred for at least one series.")
}

# ===========================
# MTM + F-test per series
# ===========================
mtm_results <- vector("list", M)

for (i in seq_len(M)) {
  x_i <- unlist(series_list$x[[i]])
  dt_i <- series_list$dt[i]
  fs_i <- series_list$fs[i]

  mtm <- multitaper::spec.mtm(ts(x_i, deltat = dt_i),
    nw = NW, k = K, Ftest = TRUE, plot = FALSE,
    deltat = dt_i, dtUnits = "s"
  )

  Fvals <- get_F(mtm)
  freqs <- mtm$freq
  pvals <- 1 - pf(Fvals, df1, df2)

  fmax_i <- if (is.null(fmax)) (fs_i / 2 - 1e-6) else min(fmax, fs_i / 2 - 1e-6)
  keep_i <- which(freqs > fmin & freqs < fmax_i)
  m_i <- length(keep_i)
  if (m_i == 0L) {
    stop(sprintf(
      "No frequency bins within (%.4g, %.4g) Hz for series %s.",
      fmin, fmax_i, as.character(series_list$PART_ID[i])
    ))
  }

  p_i <- pvals[keep_i]
  F_i <- Fvals[keep_i]
  fr_i <- freqs[keep_i]

  idx_min <- which.min(p_i)
  pmin <- p_i[idx_min]
  p_bonf <- min(1, m_i * pmin)

  mtm_results[[i]] <- list(
    PART_ID     = series_list$PART_ID[i],
    N           = length(x_i),
    dt          = dt_i,
    fs          = fs_i,
    m           = m_i,
    peak_freq   = fr_i[idx_min],
    peak_F      = F_i[idx_min],
    min_p_raw   = pmin,
    p_bonf      = p_bonf
  )
}

summary_df <- dplyr::bind_rows(lapply(mtm_results, tibble::as_tibble)) %>%
  mutate(sig_at_alpha = p_bonf <= alpha)

# ==============================
# Global Fisher combination test
# ==============================
eps <- .Machine$double.xmin
p_series <- pmax(summary_df$p_bonf, eps)
T_fisher <- -2 * sum(log(p_series))
p_global <- 1 - pchisq(T_fisher, df = 2 * M)

# ==================
# Console report
# ==================
line("=")
cat("Multitaper Harmonic F-test: Group Decision Report (Real Data)\n")
line("=")

dt_stats <- summary_df %>% summarize(
  M = dplyr::n(),
  N_min = min(N), N_median = stats::median(N), N_max = max(N),
  dt_min = min(dt), dt_median = stats::median(dt), dt_max = max(dt),
  fs_min = min(fs), fs_median = stats::median(fs), fs_max = max(fs),
  m_min = min(m), m_median = stats::median(m), m_max = max(m)
)

cat("Data settings:\n")
cat(sprintf("  Series count (M): %d\n", dt_stats$M))
cat(sprintf(
  "  Samples per series (N): min/median/max = %d / %d / %d\n",
  dt_stats$N_min, dt_stats$N_median, dt_stats$N_max
))
cat(sprintf(
  "  Sampling dt (s): min/median/max = %s / %s / %s\n",
  fmt(dt_stats$dt_min, 6), fmt(dt_stats$dt_median, 6), fmt(dt_stats$dt_max, 6)
))
cat(sprintf(
  "  Sampling fs (Hz): min/median/max = %s / %s / %s\n",
  fmt(dt_stats$fs_min, 6), fmt(dt_stats$fs_median, 6), fmt(dt_stats$fs_max, 6)
))
cat(sprintf(
  "  Frequency bins tested per series (m): min/median/max = %d / %d / %d\n",
  dt_stats$m_min, dt_stats$m_median, dt_stats$m_max
))
cat(sprintf(
  "  Frequency band scanned: f in (%.4g, %s) Hz per series\n",
  fmin, ifelse(is.null(fmax), "Nyquist", as.character(fmax))
))
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

cat("Per-series summary (Bonferroni-corrected p over frequency search):\n")
print(
  summary_df %>%
    transmute(PART_ID, N,
      dt = round(dt, 6), fs = round(fs, 6),
      peak_freq = round(peak_freq, 6),
      peak_F = round(peak_F, 3),
      min_p_raw = signif(min_p_raw, 3),
      p_bonf = signif(p_bonf, 3),
      sig_at_alpha
    ),
  row.names = FALSE
)
line("-")

# Top-K strongest series
ord <- order(summary_df$p_bonf)
kshow <- min(10, M)
cat(sprintf("Top %d series by strength (smallest per-series Bonferroni p):\n", kshow))
print(summary_df[ord[1:kshow], c("PART_ID", "peak_freq", "peak_F", "min_p_raw", "p_bonf", "sig_at_alpha")],
  row.names = FALSE
)
line("=")

# ==================
# Visualizations (APA-style)
# ==================
# Requires: summary_df, M, alpha, p_global, df1, df2, fmt() from earlier code.

library(ggplot2)

# APA-like minimalist theme
theme_apa <- function(base_size = 12, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle = element_text(size = rel(0.95)),
      plot.caption = element_text(size = rel(0.85), colour = "grey30")
    )
}

# Data for plotting
plot_df <- summary_df %>%
  mutate(
    neglog10_p = -log10(p_bonf),
    sig = ifelse(sig_at_alpha, "Significant (Bonferroni, α)", "Not significant")
  )

# (a) Histogram of per-series Bonferroni p-values with α line
p1 <- ggplot(plot_df, aes(x = p_bonf)) +
  geom_histogram(bins = 12, color = "black", fill = "grey85") +
  geom_vline(xintercept = alpha, linetype = "dashed", color = "red") +
  labs(
    title = "(a) Per-series Bonferroni-adjusted p-values",
    x = "Adjusted p-value",
    y = "Count",
    caption = paste0(
      "Red dashed line: α = ", format(alpha, digits = 2),
      "; Global Fisher p = ", fmt(p_global, 4)
    )
  ) +
  theme_apa()

# (b) QQ plot vs Uniform(0,1) with 95% Beta order-statistic bands
pp <- sort(plot_df$p_bonf)
i <- seq_along(pp)
u <- i / (M + 1)
lower <- qbeta(0.025, i, M + 1 - i)
upper <- qbeta(0.975, i, M + 1 - i)
qq_df <- tibble::tibble(u = u, p = pp, lower = lower, upper = upper)

p2 <- ggplot(qq_df, aes(x = u, y = p)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey90") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey40") +
  geom_point(shape = 16, size = 1.7) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "(b) QQ plot of per-series p-values",
    x = "Theoretical quantiles (Uniform[0,1])",
    y = "Observed adjusted p-values",
    caption = "Shaded band: pointwise 95% Beta order-statistic intervals"
  ) +
  theme_apa()

# (c) Peak frequency vs peak F with raw F critical line
Fcrit_raw <- qf(1 - alpha, df1, df2)
p3 <- ggplot(plot_df, aes(x = peak_freq, y = peak_F, shape = sig)) +
  geom_hline(yintercept = Fcrit_raw, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, color = "black") +
  scale_shape_manual(values = c("Not significant" = 1, "Significant (Bonferroni, α)" = 16)) +
  labs(
    title = "(c) Strongest line per series",
    x = "Peak frequency (Hz)",
    y = "Peak F-statistic",
    subtitle = paste0(
      "Dashed line: F critical at raw α = ", format(alpha, digits = 2),
      " (df1 = ", df1, ", df2 = ", df2, "); point fill reflects Bonferroni decision"
    )
  ) +
  theme_apa()

# (d) Top-K series by strength (−log10 adjusted p); annotate frequency and F
kshow <- min(10, M)
top_df <- plot_df %>%
  arrange(p_bonf) %>%
  dplyr::slice_head(n = kshow) %>%
  mutate(
    PART_ID = factor(PART_ID, levels = rev(PART_ID)),
    label = sprintf("%.3f Hz; F=%.2f", peak_freq, peak_F)
  )

p4 <- ggplot(top_df, aes(x = neglog10_p, y = PART_ID)) +
  geom_segment(aes(x = 0, xend = neglog10_p, y = PART_ID, yend = PART_ID), color = "grey60") +
  geom_point(aes(shape = sig), size = 2, color = "black") +
  scale_shape_manual(values = c("Not significant" = 1, "Significant (Bonferroni, α)" = 16)) +
  geom_text(aes(label = label), hjust = 0, nudge_x = 0.05, size = 3) +
  labs(
    title = "(d) Top series by strength",
    x = expression(-log[10]("adjusted p-value")),
    y = "Participant",
    caption = "Text shows peak frequency and F-statistic for the strongest line"
  ) +
  xlim(0, max(top_df$neglog10_p) * 1.15) +
  theme_apa()

# Arrange as a 2x2 grid if gridExtra is available; otherwise print sequentially
if (requireNamespace("gridExtra", quietly = TRUE)) {
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
} else {
  print(p1)
  print(p2)
  print(p3)
  print(p4)
}

# Optional: save high-resolution figures (APA-friendly)
# ggsave("fig_a_p_hist.png", p1, width = 6.5, height = 4.5, dpi = 300)
# ggsave("fig_b_qq.png",   p2, width = 6.5, height = 4.5, dpi = 300)
# ggsave("fig_c_scatter.png", p3, width = 6.5, height = 4.5, dpi = 300)
# ggsave("fig_d_topk.png", p4, width = 6.5, height = 5.0, dpi = 300)
