library(tidyverse)
library(fs)
library(patchwork)
library(lavaan)
library(ggplot2)
library(psych)
library(corrplot)
library(GPArotation)
library(semTools)
library(dplyr)      


validate_factor_analysis <- function(data, prefix = c("temporal", "inspection")) {

  df_numeric <- data |>
    select(starts_with(prefix)) |>
    drop_na()

  n_vars <- ncol(df_numeric)
  n_obs <- nrow(df_numeric)

  corr_matrix <- cor(df_numeric, use = "pairwise.complete.obs")
  kmo <- KMO(df_numeric)
  bart <- cortest.bartlett(df_numeric)
  pa <- fa.parallel(df_numeric, n.iter = 100, fa = "fa", plot = FALSE)
  efa <- fa(df_numeric, nfactors = 1, rotate = "none")

  # Safe variance explained
  var_explained <- sum(efa$Vaccounted[2]) * 100

  cat("\n=== FACTOR ANALYSIS VALIDATION ===\n")
  cat(sprintf("Variables: %d | Observations: %d\n\n", n_vars, n_obs))

  cat("CORRELATION MATRIX\n")
  print(round(corr_matrix, 2))
  cat(sprintf("Mean r = %.3f | Range = %.3f – %.3f\n\n",
              mean(corr_matrix[lower.tri(corr_matrix)]),
              min(corr_matrix[lower.tri(corr_matrix)]),
              max(corr_matrix[lower.tri(corr_matrix)])))

  cat("SAMPLING ADEQUACY\n")
  cat(sprintf("KMO = %.3f (%s)\n", kmo$MSA,
              dplyr::case_when(kmo$MSA >= 0.9 ~ "Marvelous",
                               kmo$MSA >= 0.8 ~ "Meritorious",
                               kmo$MSA >= 0.7 ~ "Middling",
                               TRUE ~ "Unacceptable")))
  cat(sprintf("Bartlett χ²(%d) = %.2f, p  = %.3f\n\n", bart$df, bart$chisq, bart$p.value))

  cat("FACTOR RETENTION\n")
  cat(sprintf("Parallel analysis suggests %d factor(s)\n\n", pa$nfact))

  cat("ONE-FACTOR SOLUTION\n")
  cat(sprintf("Eigenvalue = %.3f | Variance explained = %.1f%%\n\n",
              efa$values[1], var_explained))

  cat("LOADINGS & COMMUNALITIES\n")
  loadings_tbl <- tibble(
    Variable = names(efa$loadings[, 1]),
    Loading = round(efa$loadings[, 1], 3),
    Communality = round(efa$communality, 3)
  )
  print.data.frame(loadings_tbl, row.names = FALSE)

  invisible(list(cor = corr_matrix, kmo = kmo, bartlett = bart, pa = pa, efa = efa))
}

processing_speed <- read_csv("cleaned_cognitive_data.csv")

validate_factor_analysis(processing_speed)


## CFA

# ------------------------------------------------------------
# 1. One-factor model
# ------------------------------------------------------------
one_factor_model <- '
  processing_speed =~ temporal_SQUARES + temporal_CIRCLES +
                      inspection_CIRCLES + inspection_SQUARES
'

fit_1f <- cfa(one_factor_model,
              data = processing_speed,
              std.lv = TRUE,   # fix latent variance to 1
              estimator = "MLR")  # robust to non-normality

# ------------------------------------------------------------
# 2. Fit indices and reliability
# ------------------------------------------------------------
library(semTools)  # for reliability calculations

# Helper to pull the indices we care about
get_fit <- function(fit) {
  fitMeasures(fit, c("chisq", "df", "pvalue",
                     "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
                     "srmr"))
}

# Extract fit indices
fit1 <- get_fit(fit_1f)

# Calculate omega reliability
omega_total <- reliability(fit_1f)["omega",]
omega_hier <- reliability(fit_1f)["omega2",]  # hierarchical omega if applicable

# ------------------------------------------------------------
# 3. Enhanced Report
# ------------------------------------------------------------
cat("\n=============== ONE-FACTOR MODEL RESULTS ===============\n\n")

# Model fit
cat("MODEL FIT:\n")
cat(sprintf("χ²(%d) = %.3f, p = %.3f\n",
            fit1["df"], fit1["chisq"], fit1["pvalue"]))

# Interpretation of chi-square
if (fit1["pvalue"] > 0.05) {
  cat("✓ Non-significant χ² indicates acceptable model fit\n\n")
} else {
  cat("⚠ Significant χ² suggests model misfit (consider sample size)\n\n")
}

# Fit indices table
cat("FIT INDICES:\n")
fit_table <- data.frame(
  Index = c("CFI", "TLI", "RMSEA", "RMSEA 90% CI", "SRMR"),
  Value = c(
    sprintf("%.3f", fit1["cfi"]),
    sprintf("%.3f", fit1["tli"]),
    sprintf("%.3f", fit1["rmsea"]),
    sprintf("[%.3f, %.3f]", fit1["rmsea.ci.lower"], fit1["rmsea.ci.upper"]),
    sprintf("%.3f", fit1["srmr"])
  ),
  Cutoff = c("≥ .95", "≥ .95", "≤ .06", "Upper CI ≤ .06", "≤ .08"),
  Status = c(
    ifelse(fit1["cfi"] >= 0.95, "✓", "✗"),
    ifelse(fit1["tli"] >= 0.95, "✓", "✗"),
    ifelse(fit1["rmsea"] <= 0.06, "✓", "✗"),
    ifelse(fit1["rmsea.ci.upper"] <= 0.06, "✓", "✗"),
    ifelse(fit1["srmr"] <= 0.08, "✓", "✗")
  )
)

print(fit_table, row.names = FALSE)

# Reliability
cat("\n\nRELIABILITY:\n")
cat(sprintf("Omega Total (ωt): %.3f\n", omega_total))
if (!is.na(omega_hier)) {
  cat(sprintf("Omega Hierarchical (ωh): %.3f\n", omega_hier))
}

# Reliability interpretation
if (omega_total >= 0.80) {
  cat("✓ Excellent internal consistency (ω ≥ .80)\n")
} else if (omega_total >= 0.70) {
  cat("✓ Acceptable internal consistency (ω ≥ .70)\n")
} else {
  cat("⚠ Poor internal consistency (ω < .70)\n")
}

cat("\n===============================================\n")

# Optional: Factor loadings
cat("\nFACTOR LOADINGS:\n")
loadings_table <- parameterEstimates(fit_1f, standardized = TRUE) |>
  filter(op == "=~") |>
  select(rhs, est, std.all, pvalue) |>
  rename(Indicator = rhs,
         Unstandardized = est,
         Standardized = std.all,
         `p-value` = pvalue)

print(loadings_table, row.names = FALSE, digits = 3)






