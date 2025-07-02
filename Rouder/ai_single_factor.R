library(tidyverse)
library(fs)
library(psych)
library(lavaan)
library(corrplot)
library(semTools)

legal_ids <- c("Y2-S2-002K25", "Y2-S2-003K22", "Y2-S2-004K21", "Y2-S2-008K23",
               "Y2-S2-009K22", "Y2-S2-010K21", "Y2-S2-011K24", "Y2-S2-012M26", 
               "Y2-S2-013M31", "Y2-S2-014K21", "Y2-S2-015K27", "Y2-S2-017K23",
               "Y2-S2-018K21", "Y2-S2-019K19", "Y2-S2-020K19", "Y2-S2-021M22",
               "Y2-S2-022M27", "Y2-S2-024M20", "Y2-S2-025K21", "Y2-S2-026K24",
               "Y2-S2-027K22", "Y2-S2-028K24", "Y2-S2-029M24", "Y2-S2-030K23",
               "Y2-S2-033K24", "Y2-S2-034K26", "Y2-S2-035K29", "Y2-S2-037K20",
               "Y2-S2-045K30", "Y2-S2-047M32", "Y2-S2-049K21", "Y2-S2-051K23",
               "Y2-S2-053K19", "Y2-S2-054K22", "Y2-S2-055K22", "Y2-S2-056M23",
               "Y2-S2-058K25", "Y2-S2-062M21", "Y2-S2-065K21", "Y2-S2-068K20",
               "Y2-S2-069K23", "Y2-S2-070K25", "Y2-S2-071K29", "Y2-S2-072M23",
               "Y2-S2-074M24", "Y2-S2-076K22", "Y2-S2-077K20", "Y2-S2-078M27",
               "Y2-S2-080K20", "Y2-S2-081K20", "Y2-S2-082K20", "Y2-S2-084K38",
               "Y2-S2-088K19", "Y2-S2-090M19", "Y2-S2-092M28", "Y2-S2-094K22",
               "Y2-S2-095M27", "Y2-S2-096K22", "Y2-S2-100K21", "Y2-S2-102M22",
               "Y2-S2-104K22", "Y2-S2-107K21", "Y2-S2-108M20", "Y2-S2-109K20",
               "Y2-S2-111K20", "Y2-S2-112M21", "Y2-S2-113K19", "Y2-S2-114K19",
               "Y2-S2-118K19", "Y2-S2-119K21", "Y2-S2-120K25", "Y2-S2-121M22",
               "Y2-S2-122K22", "Y2-S2-123M22", "Y2-S2-127K20", "Y2-S2-129K31",
               "Y2-S2-130K21", "Y2-S2-131M23", "Y2-S2-132M21", "Y2-S2-133K20",
               "Y2-S2-134M23", "Y2-S2-135K23", "Y2-S2-137K20", "Y2-S2-138K24",
               "Y2-S2-140K31", "Y2-S2-141K21", "Y2-S2-142K23", "Y2-S2-144M20",
               "Y2-S2-146K21", "Y2-S2-147K19", "Y2-S2-148K21", "Y2-S2-149K18",
               "Y2-S2-150K21", "Y2-S2-152M20", "Y2-S2-153K20", "Y2-S2-154K24",
               "Y2-S2-156K20", "Y2-S2-158K20", "Y2-S2-159K18", "Y2-S2-160K20",
               "Y2-S2-161K21", "Y2-S2-164K18", "Y2-S2-166K18", "Y2-S2-169M21",
               "Y2-S2-170K23", "Y2-S2-171K21", "Y2-S2-175K21", "Y2-S2-176M26",
               "Y2-S2-179K21", "Y2-S2-181K21", "Y2-S2-183K22", "Y2-S2-184K25",
               "Y2-S2-186K19", "Y2-S2-187M32", "Y2-S2-188M25", "Y2-S2-189K21",
               "Y2-S2-190K21", "Y2-S2-191K20", "Y2-S2-193K20", "Y2-S2-194K21",
               "Y2-S2-195K21", "Y2-S2-197K23", "Y2-S2-198K23", "Y2-S2-200K18",
               "Y2-S2-202K23", "Y2-S2-203K22", "Y2-S2-206K20", "Y2-S2-207M23",
               "Y2-S2-208K19", "Y2-S2-208K26", "Y2-S2-211M23", "Y2-S2-212K24",
               "Y2-S2-213K23", "Y2-S2-214K21", "Y2-S2-215M19", "Y2-S2-216K24",
               "Y2-S2-220M19", "Y2-S2-224K20", "Y2-S2-228M20", "Y2-S2-229K21",
               "Y2-S2-230K20", "Y2-S2-231K22", "Y2-S2-232K20", "Y2-S2-233K19",
               "Y2-S2-234K20", "Y2-S2-240K20", "Y2-S2-242K22", "Y2-S2-244K18",
               "Y2-S2-247K19", "Y2-S2-250K25", "Y2-S2-251M22", "Y2-S2-253M20",
               "Y2-S2-257K25", "Y2-S2-260K19", "Y2-S2-261K24", "Y2-S2-262K26",
               "Y2-S2-265K19", "Y2-S2-266K22", "Y2-S2-267K23", "Y2-S2-268M21",
               "Y2-S2-269K20", "Y2-S2-270K21", "Y2-S2-273K18", "Y2-S2-274K23",
               "Y2-S2-276K20", "Y2-S2-278K20", "Y2-S2-280M21", "Y2-S2-281K22",
               "Y2-S2-282M23", "Y2-S2-283K20", "Y2-S2-284K20", "Y2-S2-285M22", 
               "Y2-S2-286K22", "Y2-S2-287M23")

# Create a focused single-factor analysis script
# === SINGLE FACTOR ANALYSIS FOR PROCESSING SPEED ===
# This script evaluates whether a single-factor model is appropriate

library(tidyverse)
library(psych)
library(lavaan)
library(fs)


# === 1. DATA PREPARATION ===
meaningful_revs_count <- 10

process_task_data <- function(path_pattern, prefix) {
  dir_ls(path_pattern, glob = "*.csv") |>
    map_df(read_csv, col_types = cols(
      Reversal = col_integer(),
      Reversal_count = col_integer()
    )) |>
    filter(PART_ID %in% legal_ids) |>
    filter(Reversal == 1 & 
             Reversal_count > meaningful_revs_count & 
             Training == "exp") |>
    group_by(PART_ID, Stimuli) |>
    summarize(mean_soa = mean(SOA), .groups = "drop") |>
    pivot_wider(
      names_from = Stimuli,
      values_from = mean_soa,
      names_prefix = prefix
    )
}

temporal_order <- process_task_data("./data/Temporal-order/beh", "temporal_")
inspection_time <- process_task_data("./data/Inspection-Time/beh", "inspection_")

processing_speed <- inner_join(temporal_order, inspection_time, by = "PART_ID") |> 
  drop_na()

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

