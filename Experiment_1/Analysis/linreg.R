# Packages
library(readr)
library(dplyr)
library(purrr)
library(ggplot2)
library(glmmTMB)
library(tidyverse)
library(DHARMa)
library(performance)
library(see)
library(papaja) 

# Colorblind palette
red <- "#EE6677"    # Rose/pink-red
green <- "#228833"  # Green
blue <- "#4477AA"   # Blue
purple <- "#AA3377" # Purple
cyan <- "#66CCEE"   # Cyan
yellow <- "#DDCC77" # Yellow
brown <- "#AA7744"  # Brown
gray <- "#BBBBBB"   # Gray

red <-  "#DF5069"
green <-  "#5ECF4B"
blue <-  "#3A78B8"


# 1) Load and prepare data (integrates your code; robust to multiple CSVs)
data.names <- list.files("../Data", full.names = TRUE, pattern = "*.csv$")

data.files <-
  if (length(data.names) == 1) {
    read_csv(data.names)
  } else {
    map_dfr(data.names, read_csv, .id = "source")
  }

df <- data.files |>
  filter(Trial_type == "experiment") |>
  select(PART_ID, CSI, Corr) |>
  rename(PDI = CSI, CORR = Corr) |>
  mutate(
    PDI = as.numeric(PDI) * 16.67, # Convert to ms (assuming 60 Hz refresh rate)
    PART_ID = as.factor(PART_ID),
    CORR = as.integer(CORR)
  ) |>
  filter(PDI >= 0)

## Lowess

# 2) Aggregate data across participants
df_agg <- df |>
  group_by(PDI) |>
  summarise(
    mean_CORR = mean(CORR, na.rm = TRUE),
    se_CORR = sd(CORR, na.rm = TRUE) / sqrt(n()),
    n = n()
  )

# 3) Fit LOWESS regression
lowess_fit <- lowess(df_agg$PDI, df_agg$mean_CORR, f = 0.5)

# 4) Bootstrap confidence intervals
set.seed(123)
n_boot <- 1000
boot_fits <- matrix(NA, nrow = length(lowess_fit$x), ncol = n_boot)

for (i in 1:n_boot) {
  boot_sample <- df |>
    group_by(PART_ID) |>
    slice_sample(prop = 1, replace = TRUE) |>
    ungroup() |>
    group_by(PDI) |>
    summarise(mean_CORR = mean(CORR, na.rm = TRUE))
  
  boot_lowess <- lowess(boot_sample$PDI, boot_sample$mean_CORR, f = 0.5)
  boot_fits[, i] <- approx(boot_lowess$x, boot_lowess$y, xout = lowess_fit$x)$y
}

# Calculate confidence intervals
ci_lower <- apply(boot_fits, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_upper <- apply(boot_fits, 1, quantile, probs = 0.975, na.rm = TRUE)

# 5) Create plot
lowess_df <- tibble(
  PDI = lowess_fit$x,
  fitted = lowess_fit$y,
  ci_lower = ci_lower,
  ci_upper = ci_upper
)

p1 <- ggplot() +
  # geom_point(data = df_agg, aes(x = PDI, y = mean_CORR), alpha = 0.5) +
  geom_ribbon(
    data = lowess_df,
    aes(x = PDI, ymin = ci_lower, ymax = ci_upper),
    alpha = 0.3, fill = blue
  ) +
  geom_smooth(
    data = df_agg,
    aes(x = PDI, y = mean_CORR),
    method = "lm", color = red,
    se = FALSE, linetype = "dashed"
  ) +
  geom_line(
    data = lowess_df,
    aes(x = PDI, y = fitted),
    color = blue, linewidth = 1
  ) +
  labs(
    x = "PDI (ms)",
    y = "Mean Antisaccade Accuracy",
    title = "LOWESS Regression: PDI vs Accuracy"
  ) +
  papaja::theme_apa() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )


ggsave("lowess.png",p1,  width = 10, height = 6, dpi = 600)

## 



# 2) Aggregate to counts per participant x PDI
df_cell <- df |>
  group_by(PART_ID, PDI) |>
  summarise(
    correct_n = sum(CORR, na.rm = TRUE),
    total_n = n(),
    .groups = "drop"
  )

# Sanity check: should be 11 trials per cell
if (any(df_cell$total_n != 11)) {
  warning("Some PART_ID × PDI cells do not have 11 trials. Using actual total_n in the model.")
}

# 3) Predictors: linear (per 100 ms), quadratic, cubic, and log (log1p for PDI=0 safety)
df_cell <- df_cell |>
  mutate(
    PDI100  = (PDI - mean(PDI, na.rm = TRUE)) / 100,
    PDI1002 = PDI100^2,
    PDI1003 = PDI100^3,
    lPDIc   = log1p(PDI) - mean(log1p(PDI), na.rm = TRUE) # centered log scale
  )

# 4) Models (keep random slope for the main linear/log predictor; higher-order terms fixed only)
# Linear
fit_lin <- glmmTMB(
  cbind(correct_n, total_n - correct_n) ~ PDI100 + (1 + PDI100 | PART_ID),
  family = binomial(link = "logit"),
  data = df_cell
)
if (!is.null(fit_lin$sdr$pdHess) && !fit_lin$sdr$pdHess) {
  message("Linear: refitting with random intercept only...")
  fit_lin <- glmmTMB(
    cbind(correct_n, total_n - correct_n) ~ PDI100 + (1 | PART_ID),
    family = binomial(link = "logit"),
    data = df_cell
  )
}

# Quadratic
fit_quad <- glmmTMB(
  cbind(correct_n, total_n - correct_n) ~ PDI100 + PDI1002 + (1 + PDI100 | PART_ID),
  family = binomial(link = "logit"),
  data = df_cell
)
if (!is.null(fit_quad$sdr$pdHess) && !fit_quad$sdr$pdHess) {
  message("Quadratic: refitting with random intercept only...")
  fit_quad <- glmmTMB(
    cbind(correct_n, total_n - correct_n) ~ PDI100 + PDI1002 + (1 | PART_ID),
    family = binomial(link = "logit"),
    data = df_cell
  )
}

# Cubic
fit_cubic <- glmmTMB(
  cbind(correct_n, total_n - correct_n) ~ PDI100 + PDI1002 + PDI1003 + (1 + PDI100 | PART_ID),
  family = binomial(link = "logit"),
  data = df_cell
)
if (!is.null(fit_cubic$sdr$pdHess) && !fit_cubic$sdr$pdHess) {
  message("Cubic: refitting with random intercept only...")
  fit_cubic <- glmmTMB(
    cbind(correct_n, total_n - correct_n) ~ PDI100 + PDI1002 + PDI1003 + (1 | PART_ID),
    family = binomial(link = "logit"),
    data = df_cell
  )
}

# Log
fit_log <- glmmTMB(
  cbind(correct_n, total_n - correct_n) ~ lPDIc + (1 + lPDIc | PART_ID),
  family = binomial(link = "logit"),
  data = df_cell
)
if (!is.null(fit_log$sdr$pdHess) && !fit_log$sdr$pdHess) {
  message("Log: refitting with random intercept only...")
  fit_log <- glmmTMB(
    cbind(correct_n, total_n - correct_n) ~ lPDIc + (1 | PART_ID),
    family = binomial(link = "logit"),
    data = df_cell
  )
}

# 5) Formal tests
## (a) Linear slope vs. no-trend
fit0 <- update(fit_lin, . ~ . - PDI100)
lrt_lin <- as.data.frame(anova(fit0, fit_lin))
cat("\n=== Linear Trend Test (population-level) ===\n")
coefs_lin <- fixef(fit_lin)$cond
vc_lin <- vcov(fit_lin)$cond
b_lin <- unname(coefs_lin["PDI100"])
se_lin <- sqrt(vc_lin["PDI100", "PDI100"])
z_lin <- b_lin / se_lin
p_lin <- 2 * pnorm(abs(z_lin), lower.tail = FALSE)
OR_lin <- exp(b_lin)
OR_lo <- exp(b_lin - 1.96 * se_lin)
OR_hi <- exp(b_lin + 1.96 * se_lin)
cat(sprintf("Log-odds slope per +100 ms:  b = %.4f, SE = %.4f, z = %.2f, p = %.4g\n", b_lin, se_lin, z_lin, p_lin))
cat(sprintf("Odds ratio per +100 ms:      OR = %.3f, 95%% CI [%.3f, %.3f]\n", OR_lin, OR_lo, OR_hi))
cat(sprintf(
  "LRT vs. no-trend:            Chi^2 = %.2f on %d df, p = %.4g\n\n",
  lrt_lin$Chisq[2], lrt_lin$Df[2], lrt_lin$`Pr(>Chisq)`[2]
))

## (b) Within-polynomial family (nested): Linear vs Quadratic vs Cubic
lrt_lin_quad <- as.data.frame(anova(fit_lin, fit_quad))
lrt_quad_cub <- as.data.frame(anova(fit_quad, fit_cubic))

coefs_quad <- fixef(fit_quad)$cond
b_quad2 <- unname(coefs_quad["PDI1002"])
se_quad2 <- sqrt(vcov(fit_quad)$cond["PDI1002", "PDI1002"])
z_quad2 <- b_quad2 / se_quad2
p_quad2 <- 2 * pnorm(abs(z_quad2), lower.tail = FALSE)

coefs_cub <- fixef(fit_cubic)$cond
b_cub3 <- unname(coefs_cub["PDI1003"])
se_cub3 <- sqrt(vcov(fit_cubic)$cond["PDI1003", "PDI1003"])
z_cub3 <- b_cub3 / se_cub3
p_cub3 <- 2 * pnorm(abs(z_cub3), lower.tail = FALSE)

cat("=== Polynomial Family: Nested LRTs ===\n")
cat(sprintf(
  "Quadratic vs Linear:         Chi^2 = %.2f on %d df, p = %.4g\n",
  lrt_lin_quad$Chisq[2], lrt_lin_quad$Df[2], lrt_lin_quad$`Pr(>Chisq)`[2]
))
cat(sprintf(
  "Cubic vs Quadratic:          Chi^2 = %.2f on %d df, p = %.4g\n\n",
  lrt_quad_cub$Chisq[2], lrt_quad_cub$Df[2], lrt_quad_cub$`Pr(>Chisq)`[2]
))
cat(sprintf(
  "Quadratic term z-test:       b2 = %.4f, SE = %.4f, z = %.2f, p = %.4g\n",
  b_quad2, se_quad2, z_quad2, p_quad2
))
cat(sprintf(
  "Cubic term z-test:           b3 = %.4f, SE = %.4f, z = %.2f, p = %.4g\n\n",
  b_cub3, se_cub3, z_cub3, p_cub3
))

## (c) Non-nested comparison including log model: use AIC
aic_tab <- AIC(fit_lin, fit_quad, fit_cubic, fit_log)
aic_tab <- aic_tab[order(aic_tab$AIC), ]
deltaAIC <- aic_tab$AIC - min(aic_tab$AIC)
cat("=== AIC Comparison (non-nested) ===\n")
for (i in seq_len(nrow(aic_tab))) {
  cat(sprintf(
    "%-12s AIC = %.1f, ΔAIC = %.1f\n",
    rownames(aic_tab)[i], aic_tab$AIC[i], deltaAIC[i]
  ))
}
cat(sprintf("\n=> Best by AIC: %s\n\n", rownames(aic_tab)[1]))

# 6) Visualization: all four population-level fits with 95% CI
grid <- data.frame(PDI = seq(min(df_cell$PDI), max(df_cell$PDI), length.out = 300)) |>
  mutate(
    PDI100  = (PDI - mean(df_cell$PDI, na.rm = TRUE)) / 100,
    PDI1002 = PDI100^2,
    PDI1003 = PDI100^3,
    lPDIc   = log1p(PDI) - mean(log1p(df_cell$PDI), na.rm = TRUE)
  )

# Predictions on link scale, no random effects (population-average)
pred_lin <- predict(fit_lin, newdata = grid, type = "link", se.fit = TRUE, re.form = NA)
pred_quad <- predict(fit_quad, newdata = grid, type = "link", se.fit = TRUE, re.form = NA)
pred_cubic <- predict(fit_cubic, newdata = grid, type = "link", se.fit = TRUE, re.form = NA)
pred_log <- predict(fit_log, newdata = grid, type = "link", se.fit = TRUE, re.form = NA)

make_df <- function(PDI, pr, label) {
  data.frame(
    PDI = PDI,
    fit = plogis(pr$fit),
    lo = plogis(pr$fit - 1.96 * pr$se.fit),
    hi = plogis(pr$fit + 1.96 * pr$se.fit),
    model = label
  )
}

grid_long <- bind_rows(
  make_df(grid$PDI, pred_lin, "Linear"),
  make_df(grid$PDI, pred_quad, "Quadratic"),
  make_df(grid$PDI, pred_cubic, "Cubic"),
  make_df(grid$PDI, pred_log, "Log")
)

p_all <- ggplot(grid_long, aes(x = PDI, y = fit, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Linear" = "#2c7fb8", "Quadratic" = "#f03b20", "Cubic" = "#6a3d9a", "Log" = "#33a02c")) +
  scale_fill_manual(values = c("Linear" = "#2c7fb8", "Quadratic" = "#f03b20", "Cubic" = "#6a3d9a", "Log" = "#33a02c")) +
  labs(
    x = "PDI (ms)",
    y = "Accuracy",
    title = "Accuracy vs. PDI: Model fits (population-average with 95% CI)"
  ) +
  theme_classic()

print(p_all)


library(DHARMa)

set.seed(1)
sim <- simulateResiduals(fit_log, n = 1000, plot = FALSE)

# One multi-panel diagnostic figure
plot(sim)

# Safe p-value extractor
get_p <- function(x) {
  if (is.null(x)) return(NA_real_)
  p <- tryCatch(x$p.value, error = function(e) NA_real_)
  if (is.null(p)) return(NA_real_)
  as.numeric(p)[1]
}

# Key tests (robust to missing/NA)
omni <- tryCatch(testResiduals(sim, plot = FALSE),     error = function(e) NULL)
uni  <- tryCatch(testUniformity(sim, plot = FALSE),    error = function(e) NULL)
disp <- tryCatch(testDispersion(sim, plot = FALSE),    error = function(e) NULL)
zi   <- tryCatch(testZeroInflation(sim, plot = FALSE), error = function(e) NULL)
out  <- tryCatch(testOutliers(sim, plot = FALSE, type = "bootstrap"),      error = function(e) NULL)

res_tests <- data.frame(
  test    = c("Omnibus", "Uniformity", "Dispersion", "Zero-inflation", "Outliers"),
  p_value = signif(c(get_p(omni), get_p(uni), get_p(disp), get_p(zi), get_p(out)), 3)
)

print(res_tests, row.names = FALSE)