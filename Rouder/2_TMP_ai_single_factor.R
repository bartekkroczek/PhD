library(tidyverse)
library(fs)
library(psych)
library(lavaan)
library(corrplot)
library(semTools)

processing_speed <- read_csv("cleaned_cognitive_data.csv")

df_numeric <- processing_speed |> 
  select(temporal_CIRCLES, temporal_SQUARES, inspection_CIRCLES, inspection_SQUARES)

cat("Sample size:", nrow(processing_speed), "\n\n")

# === 2. SAMPLING ADEQUACY TESTS ===
cat("=== SAMPLING ADEQUACY FOR FACTOR ANALYSIS ===\n\n")

# KMO Test
kmo_result <- KMO(df_numeric)
cat("Kaiser-Meyer-Olkin (KMO) Measure:", round(kmo_result$MSA, 3), "\n")

if (kmo_result$MSA < 0.6) {
  cat("⚠️ WARNING: KMO < 0.6 suggests factor analysis may be inappropriate\n")
} else {
  cat("✓ KMO indicates adequate sampling for factor analysis\n")
}

# Bartlett's Test
bart_result <- cortest.bartlett(df_numeric)
cat("\nBartlett's Test p-value:", format.pval(bart_result$p.value, digits = 3), "\n")

if (bart_result$p.value < 0.05) {
  cat("✓ Variables are sufficiently correlated for factor analysis\n")
} else {
  cat("⚠️ WARNING: Variables may be too uncorrelated for factor analysis\n")
}

# === 3. EXPLORATORY FACTOR ANALYSIS - SINGLE FACTOR ===
cat("\n=== SINGLE FACTOR EXPLORATORY ANALYSIS ===\n")

# Fit single factor model
efa_1f <- fa(df_numeric, nfactors = 1, rotate = "none", fm = "ml")

# Display loadings
cat("\nFactor Loadings:\n")
print(round(efa_1f$loadings, 3))

# Variance explained
var_explained <- sum(efa_1f$communality) / length(efa_1f$communality) * 100
cat("\nVariance Explained:", round(var_explained, 1), "%\n")

# Fit statistics
cat("\nModel Fit (EFA):\n")
cat("- Chi-square:", round(efa_1f$STATISTIC, 2), "\n")
cat("- df:", efa_1f$dof, "\n")
cat("- p-value:", format.pval(efa_1f$PVAL, digits = 3), "\n")
cat("- RMSEA:", round(efa_1f$RMSEA[1], 3), "\n")
cat("- TLI:", round(efa_1f$TLI, 3), "\n")

# === 4. CONFIRMATORY FACTOR ANALYSIS - SINGLE FACTOR ===
cat("\n=== SINGLE FACTOR CONFIRMATORY ANALYSIS ===\n")

model_1factor <- '
  processing_speed =~ temporal_CIRCLES + temporal_SQUARES + 
                      inspection_CIRCLES + inspection_SQUARES
'

fit_1f <- cfa(model_1factor, data = processing_speed, std.lv = TRUE)

# Extract fit measures
fit_measures <- fitMeasures(fit_1f, c("chisq", "df", "pvalue", "cfi", "tli", 
                                      "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
                                      "srmr"))

cat("\nModel Fit (CFA):\n")
cat("- Chi-square:", round(fit_measures["chisq"], 2), "\n")
cat("- df:", fit_measures["df"], "\n")
cat("- p-value:", format.pval(fit_measures["pvalue"], digits = 3), "\n")
cat("- CFI:", round(fit_measures["cfi"], 3), "\n")
cat("- TLI:", round(fit_measures["tli"], 3), "\n")
cat("- RMSEA:", round(fit_measures["rmsea"], 3), 
    "[", round(fit_measures["rmsea.ci.lower"], 3), ",", 
    round(fit_measures["rmsea.ci.upper"], 3), "]\n")
cat("- SRMR:", round(fit_measures["srmr"], 3), "\n")

# Standardized loadings
std_solution <- standardizedSolution(fit_1f)
loadings <- std_solution[std_solution$op == "=~", c("rhs", "est.std")]
cat("\nStandardized Loadings (CFA):\n")
print(loadings)

# === 5. RELIABILITY ANALYSIS ===
cat("\n=== RELIABILITY OF SINGLE FACTOR ===\n")

# Cronbach's alpha
alpha_result <- alpha(df_numeric)
cat("Cronbach's Alpha:", round(alpha_result$total$raw_alpha, 3), "\n")

# Average inter-item correlation
avg_r <- alpha_result$total$average_r
cat("Average inter-item correlation:", round(avg_r, 3), "\n")

# Item-total correlations
item_total <- alpha_result$item.stats[, "r.drop"]
cat("\nItem-total correlations:\n")
print(round(item_total, 3))

# === 6. RESIDUAL ANALYSIS ===
cat("\n=== RESIDUAL CORRELATIONS ===\n")

# Get residual correlations - correct method
resid_output <- lavResiduals(fit_1f)
residuals <- resid_output$cov  # This is the residual covariance/correlation matrix

cat("Residual correlation matrix:\n")
print(round(residuals, 3))

# Check for large residuals
residuals_abs <- abs(residuals)
diag(residuals_abs) <- 0
large_resid <- which(residuals_abs > 0.1, arr.ind = TRUE)

if (nrow(large_resid) > 0) {
  cat("\n⚠️ Large residual correlations (|r| > 0.1):\n")
  for (i in 1:nrow(large_resid)) {
    if (large_resid[i, 1] < large_resid[i, 2]) {
      cat(sprintf("  %s - %s: %.3f\n",
                  rownames(residuals)[large_resid[i, 1]],
                  colnames(residuals)[large_resid[i, 2]],
                  residuals[large_resid[i, 1], large_resid[i, 2]]))
    }
  }
} else {
  cat("\n✓ No large residual correlations found\n")
}

# === 7. EVALUATION CRITERIA ===
cat("\n=== SINGLE FACTOR MODEL EVALUATION ===\n")
cat("=====================================\n")

# Define evaluation criteria
good_fit <- list(
  cfi = fit_measures["cfi"] > 0.95,
  tli = fit_measures["tli"] > 0.95,
  rmsea = fit_measures["rmsea"] < 0.06,
  srmr = fit_measures["srmr"] < 0.08,
  alpha = alpha_result$total$raw_alpha > 0.70,
  variance = var_explained > 50
)

# Summary evaluation
cat("\nFit Criteria Assessment:\n")
cat("- CFI > 0.95:", ifelse(good_fit$cfi, "✓ Pass", "✗ Fail"), 
    "(", round(fit_measures["cfi"], 3), ")\n")
cat("- TLI > 0.95:", ifelse(good_fit$tli, "✓ Pass", "✗ Fail"), 
    "(", round(fit_measures["tli"], 3), ")\n")
cat("- RMSEA < 0.06:", ifelse(good_fit$rmsea, "✓ Pass", "✗ Fail"), 
    "(", round(fit_measures["rmsea"], 3), ")\n")
cat("- SRMR < 0.08:", ifelse(good_fit$srmr, "✓ Pass", "✗ Fail"), 
    "(", round(fit_measures["srmr"], 3), ")\n")
cat("- Alpha > 0.70:", ifelse(good_fit$alpha, "✓ Pass", "✗ Fail"), 
    "(", round(alpha_result$total$raw_alpha, 3), ")\n")
cat("- Variance > 50%:", ifelse(good_fit$variance, "✓ Pass", "✗ Fail"), 
    "(", round(var_explained, 1), "%)\n")

# === 8. FINAL RECOMMENDATION ===
cat("\n=== CONCLUSION: IS SINGLE FACTOR APPROPRIATE? ===\n")
cat("================================================\n")

n_criteria_met <- sum(unlist(good_fit))

if (n_criteria_met >= 5) {
  cat("\n✅ RECOMMENDATION: Single factor model is APPROPRIATE\n")
  cat("- Most fit criteria are met\n")
  cat("- Items show adequate internal consistency\n")
  cat("- Single processing speed score can be used\n")
} else if (n_criteria_met >= 3) {
  cat("\n⚠️ RECOMMENDATION: Single factor model is MARGINAL\n")
  cat("- Some fit criteria are not met\n")
  cat("- Consider alternative models or item-level analysis\n")
  cat("- Use single score with caution\n")
} else {
  cat("\n❌ RECOMMENDATION: Single factor model is NOT APPROPRIATE\n")
  cat("- Most fit criteria are failed\n")
  cat("- Items may measure different constructs\n")
  cat("- Consider multi-factor models or separate analyses\n")
}

# Specific issues
cat("\nSpecific Concerns:\n")
if (!good_fit$alpha) {
  cat("- Low reliability (α =", round(alpha_result$total$raw_alpha, 3), 
      ") indicates high measurement error\n")
}
if (!good_fit$variance) {
  cat("- Low variance explained (", round(var_explained, 1), 
      "%) suggests missing factors\n")
}
if (nrow(large_resid) > 0) {
  cat("- Large residual correlations indicate local dependencies\n")
}
if (avg_r < 0.3) {
  cat("- Low average inter-item correlation (r =", round(avg_r, 3), 
      ") questions unidimensionality\n")
}

# Alternative recommendations
cat("\nAlternative Approaches:\n")
if (!good_fit$alpha || avg_r < 0.3) {
  cat("→ Analyze items separately rather than creating composite\n")
  cat("→ Use best single indicator per task\n")
  cat("→ Consider measurement model with correlated errors\n")
}

