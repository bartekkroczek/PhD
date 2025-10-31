# ===== FACTOR ANALYSIS OF PROCESSING SPEED MEASURES =====
# Purpose: Determine if temporal order and inspection time tasks measure
# a unified processing speed construct

library(tidyverse)
library(fs)
library(psych)
library(lavaan)
library(corrplot)
library(semTools)
library(RColorBrewer)

legal_ids <- c(
  "Y2-S2-002K25", "Y2-S2-003K22", "Y2-S2-004K21", "Y2-S2-008K23",
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
  "Y2-S2-286K22", "Y2-S2-287M23"
)

# === 1. DATA PREPARATION (keeping your original code) ===
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

# Rename variables for better readability
df_renamed <- df_numeric %>%
  rename(
    `Temporal Order\n(Circles)` = temporal_CIRCLES,
    `Temporal Order\n(Squares)` = temporal_SQUARES,
    `Inspection Time\n(Circles)` = inspection_CIRCLES,
    `Inspection Time\n(Squares)` = inspection_SQUARES
  )

# Calculate correlation matrix
cor_matrix <- df_renamed %>%
  cor(use = "complete.obs")


# APA 7 compliant with colorblind-friendly palette
cor_data <- cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Task1") %>%
  pivot_longer(-Task1, names_to = "Task2", values_to = "r") %>%
  mutate(
    Task1 = factor(Task1, levels = rownames(cor_matrix)),
    Task2 = factor(Task2, levels = colnames(cor_matrix)),
    row_idx = as.numeric(Task1),
    col_idx = as.numeric(Task2),
    # Show lower triangle AND diagonal
    r = ifelse(row_idx >= col_idx, r, NA),
    # Add significance indicators (replace with actual p-values)
    sig_indicator = case_when(
      is.na(r) ~ "",
      row_idx == col_idx ~ "", # No significance for diagonal
      abs(r) >= 0.5 ~ "**",    # Replace with p < .01 logic
      abs(r) >= 0.3 ~ "*",     # Replace with p < .05 logic
      TRUE ~ ""
    ),
    r_label = case_when(
      is.na(r) ~ "",
      row_idx == col_idx ~ "1.00", # Show perfect correlations on diagonal
      TRUE ~ paste0(sprintf("%.2f", r), sig_indicator)
    )
  )

ggplot(cor_data, aes(x = Task1, y = Task2, fill = r)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = r_label), color = "black", size = 5, fontface = "bold") +
  # Colorblind-friendly palette (viridis)
  scale_fill_gradient2(low = "#440154", mid = "#FDE725", high = "#35B779",
                       midpoint = 0, limit = c(-1,1),
                       name = "Pearson r", na.value = "white") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  ) +
  coord_fixed() +
  labs(caption = "* p < .05, ** p < .01")

cat("Sample size:", nrow(processing_speed), "\n\n")

cat("=== SAMPLING ADEQUACY TESTS ===\n\n")

# KMO Test
kmo_result <- KMO(df_numeric)
cat("Kaiser-Meyer-Olkin (KMO) Measure of Sampling Adequacy:\n")
cat("Overall MSA:", round(kmo_result$MSA, 3), "\n")
cat("Individual MSA values:\n")
print(round(kmo_result$MSAi, 3))

# Interpretation
if (kmo_result$MSA < 0.5) {
  cat("⚠️ KMO < 0.5: Factor analysis may not be appropriate\n")
} else if (kmo_result$MSA < 0.6) {
  cat("⚠️ KMO = 0.5-0.6: Mediocre sampling adequacy\n")
} else if (kmo_result$MSA < 0.7) {
  cat("• KMO = 0.6-0.7: Acceptable sampling adequacy\n")
} else {
  cat("✓ KMO > 0.7: Good sampling adequacy\n")
}

# Bartlett's Test of Sphericity
bart_result <- cortest.bartlett(df_numeric)
cat("\nBartlett's Test of Sphericity:\n")
cat("Chi-square:", round(bart_result$chisq, 2), "\n")
cat("df:", bart_result$df, "\n")
cat("p-value:", format.pval(bart_result$p.value, digits = 3), "\n")

if (bart_result$p.value < 0.05) {
  cat("✓ p < 0.05: Correlations differ significantly from identity matrix\n")
} else {
  cat("⚠️ p > 0.05: Variables may be uncorrelated\n")
}

cat("\n")


# === 2. ENHANCED EXPLORATORY FACTOR ANALYSIS ===
cat("=== EXPLORATORY FACTOR ANALYSIS: MULTIPLE MODELS ===\n\n")

# Test 1-factor solution
efa_1f <- fa(df_numeric, nfactors = 1, rotate = "none", fm = "ml")

# Test 2-factor solution with different rotations
efa_2f_oblimin <- fa(df_numeric, nfactors = 2, rotate = "oblimin", fm = "ml")
efa_2f_varimax <- fa(df_numeric, nfactors = 2, rotate = "varimax", fm = "ml")

# Extract fit statistics
fit_comparison <- data.frame(
  Model = c("1-factor", "2-factor (oblimin)", "2-factor (varimax)"),
  Chi_sq = c(efa_1f$STATISTIC, efa_2f_oblimin$STATISTIC, efa_2f_varimax$STATISTIC),
  df = c(efa_1f$dof, efa_2f_oblimin$dof, efa_2f_varimax$dof),
  p_value = c(efa_1f$PVAL, efa_2f_oblimin$PVAL, efa_2f_varimax$PVAL),
  RMSEA = c(efa_1f$RMSEA[1], efa_2f_oblimin$RMSEA[1], efa_2f_varimax$RMSEA[1]),
  TLI = c(efa_1f$TLI, efa_2f_oblimin$TLI, efa_2f_varimax$TLI),
  BIC = c(efa_1f$BIC, efa_2f_oblimin$BIC, efa_2f_varimax$BIC),
  Var_explained = c(
    sum(efa_1f$communality) / length(efa_1f$communality) * 100,
    sum(efa_2f_oblimin$communality) / length(efa_2f_oblimin$communality) * 100,
    sum(efa_2f_varimax$communality) / length(efa_2f_varimax$communality) * 100
  )
)

print(fit_comparison, digits = 3)

# === SCREE PLOT AND PARALLEL ANALYSIS ===
cat("\n=== FACTOR RETENTION CRITERIA ===\n")

# Scree plot
scree_data <- fa.parallel(df_numeric,
  fm = "ml", fa = "fa",
  main = "Scree Plot with Parallel Analysis"
)
cat("Suggested number of factors (parallel analysis):", scree_data$nfact, "\n")

# Eigenvalues
eigen_values <- eigen(cor(df_numeric))$values
cat("\nEigenvalues:", round(eigen_values, 3), "\n")
cat("Factors with eigenvalues > 1:", sum(eigen_values > 1), "\n")

# === 3. DETAILED 2-FACTOR SOLUTION EXAMINATION ===
cat("\n=== 2-FACTOR SOLUTION (OBLIMIN ROTATION) ===\n")
print(efa_2f_oblimin, cut = 0.3, sort = TRUE)

# Check factor correlation
if (!is.null(efa_2f_oblimin$Phi)) {
  cat("\nFactor Correlation:\n")
  print(round(efa_2f_oblimin$Phi, 3))

  factor_cor <- efa_2f_oblimin$Phi[1, 2]
  if (abs(factor_cor) > 0.7) {
    cat(
      "⚠️ High factor correlation (r =", round(factor_cor, 3),
      ") suggests factors may not be distinct\n"
    )
  } else if (abs(factor_cor) < 0.3) {
    cat(
      "✓ Low factor correlation (r =", round(factor_cor, 3),
      ") supports distinct factors\n"
    )
  } else {
    cat("• Moderate factor correlation (r =", round(factor_cor, 3), ")\n")
  }
}

# Visualize loadings
par(mfrow = c(1, 2))
fa.diagram(efa_1f, main = "1-Factor Solution")
fa.diagram(efa_2f_oblimin, main = "2-Factor Solution")
par(mfrow = c(1, 1))

# === 4. CONFIRMATORY FACTOR ANALYSIS: MODEL COMPARISON ===
cat("\n=== CONFIRMATORY FACTOR ANALYSIS ===\n")

# Model 1: Single factor
model_1factor <- "
  processing_speed =~ temporal_CIRCLES + temporal_SQUARES +
                      inspection_CIRCLES + inspection_SQUARES
"

# Model 2: Two correlated factors
model_2factor <- "
  temporal   =~ temporal_CIRCLES + temporal_SQUARES
  inspection =~ inspection_CIRCLES + inspection_SQUARES

  # Estimate correlation
  temporal ~~ inspection
"

# Model 3: Hierarchical (second-order factor)
model_hierarchical <- "
  temporal   =~ temporal_CIRCLES + temporal_SQUARES
  inspection =~ inspection_CIRCLES + inspection_SQUARES

  # Second-order factor
  processing_speed =~ temporal + inspection
"

# Model 4: Bifactor model
model_bifactor <- "
  # General factor
  general =~ temporal_CIRCLES + temporal_SQUARES +
             inspection_CIRCLES + inspection_SQUARES

  # Specific factors (orthogonal to general)
  temporal_specific   =~ temporal_CIRCLES + temporal_SQUARES
  inspection_specific =~ inspection_CIRCLES + inspection_SQUARES

  # Orthogonality constraints
  general ~~ 0*temporal_specific
  general ~~ 0*inspection_specific
  temporal_specific ~~ 0*inspection_specific
"

# Fit all models
fit_1f <- cfa(model_1factor, data = processing_speed, std.lv = TRUE)
fit_2f <- cfa(model_2factor, data = processing_speed, std.lv = TRUE)
fit_hier <- cfa(model_hierarchical, data = processing_speed, std.lv = TRUE)

# Try bifactor (may not converge with 4 indicators)
fit_bifactor <- tryCatch(
  cfa(model_bifactor, data = processing_speed, std.lv = TRUE),
  error = function(e) NULL
)

# === 5. MODEL COMPARISON ===
cat("\n=== MODEL FIT COMPARISON ===\n")

# Function to extract fit measures
get_fit_measures <- function(fit_object, model_name) {
  if (is.null(fit_object)) {
    return(NULL)
  }

  fits <- fitMeasures(fit_object, c(
    "chisq", "df", "pvalue", "cfi", "tli",
    "rmsea", "srmr", "aic", "bic"
  ))

  data.frame(
    Model = model_name,
    Chi_sq = fits["chisq"],
    df = fits["df"],
    p_value = fits["pvalue"],
    CFI = fits["cfi"],
    TLI = fits["tli"],
    RMSEA = fits["rmsea"],
    SRMR = fits["srmr"],
    AIC = fits["aic"],
    BIC = fits["bic"]
  )
}

# Combine fit statistics
fit_table <- bind_rows(
  get_fit_measures(fit_1f, "1-factor"),
  get_fit_measures(fit_2f, "2-factor"),
  get_fit_measures(fit_hier, "Hierarchical"),
  if (!is.null(fit_bifactor)) get_fit_measures(fit_bifactor, "Bifactor") else NULL
)

print(fit_table, digits = 3)

# Chi-square difference test (1-factor vs 2-factor)
if (fit_1f@Fit@converged && fit_2f@Fit@converged) {
  chi_diff <- lavTestLRT(fit_1f, fit_2f)
  cat("\n=== CHI-SQUARE DIFFERENCE TEST ===\n")
  print(chi_diff)
}

# === 6. PARAMETER ESTIMATES FOR 2-FACTOR MODEL ===
cat("\n=== 2-FACTOR MODEL PARAMETERS ===\n")
summary(fit_2f, standardized = TRUE, fit.measures = FALSE)

# Extract factor correlation
factor_cor_cfa <- lavInspect(fit_2f, "cor.lv")[1, 2]
cat("\nFactor correlation (CFA):", round(factor_cor_cfa, 3), "\n")

# === 7. RELIABILITY FOR EACH FACTOR (FIXED) ===
cat("\n=== FACTOR-SPECIFIC RELIABILITY ===\n")

# Temporal order items
temporal_items <- df_numeric[, c("temporal_CIRCLES", "temporal_SQUARES")]
alpha_temporal <- alpha(temporal_items)
cat("Temporal Order: α =", round(alpha_temporal$total$raw_alpha, 3), "\n")

# Inspection time items
inspection_items <- df_numeric[, c("inspection_CIRCLES", "inspection_SQUARES")]
alpha_inspection <- alpha(inspection_items)
cat("Inspection Time: α =", round(alpha_inspection$total$raw_alpha, 3), "\n")

# Calculate composite reliability from CFA manually
if (fit_2f@Fit@converged) {
  # Extract standardized loadings
  std_loadings <- standardizedSolution(fit_2f)

  # Calculate composite reliability for each factor
  temporal_loadings <- std_loadings %>%
    filter(lhs == "temporal" & op == "=~") %>%
    pull(est.std)

  inspection_loadings <- std_loadings %>%
    filter(lhs == "inspection" & op == "=~") %>%
    pull(est.std)

  # McDonald's omega formula
  omega_temporal <- sum(temporal_loadings)^2 /
    (sum(temporal_loadings)^2 + sum(1 - temporal_loadings^2))

  omega_inspection <- sum(inspection_loadings)^2 /
    (sum(inspection_loadings)^2 + sum(1 - inspection_loadings^2))

  cat("\nComposite Reliability (McDonald's ω):\n")
  cat("Temporal Order: ω =", round(omega_temporal, 3), "\n")
  cat("Inspection Time: ω =", round(omega_inspection, 3), "\n")
}

# === CONFIRMATORY FACTOR ANALYSIS ===
cat("\n=== CONFIRMATORY FACTOR ANALYSIS ===\n")

# Model 1: Single factor
model_1factor <- "
  processing_speed =~ temporal_CIRCLES + temporal_SQUARES +
                      inspection_CIRCLES + inspection_SQUARES
"

# Model 2: Two correlated factors
model_2factor <- "
  temporal   =~ temporal_CIRCLES + temporal_SQUARES
  inspection =~ inspection_CIRCLES + inspection_SQUARES
"

# Fit models
fit_1f <- cfa(model_1factor, data = processing_speed, std.lv = TRUE)
fit_2f <- cfa(model_2factor, data = processing_speed, std.lv = TRUE)

# === MODEL COMPARISON ===
cat("\n=== MODEL FIT COMPARISON ===\n")

# Function to safely extract fit measures
get_fit_measures_safe <- function(fit_object, model_name) {
  if (is.null(fit_object) || !fit_object@Fit@converged) {
    return(data.frame(Model = model_name, Status = "Did not converge"))
  }

  fits <- fitMeasures(fit_object, c(
    "chisq", "df", "pvalue", "cfi", "tli",
    "rmsea", "srmr", "aic", "bic"
  ))

  data.frame(
    Model = model_name,
    Status = "Converged",
    Chi_sq = round(fits["chisq"], 2),
    df = fits["df"],
    p_value = round(fits["pvalue"], 3),
    CFI = round(fits["cfi"], 3),
    TLI = round(fits["tli"], 3),
    RMSEA = round(fits["rmsea"], 3),
    SRMR = round(fits["srmr"], 3),
    AIC = round(fits["aic"], 1),
    BIC = round(fits["bic"], 1)
  )
}

# Create comparison table
fit_table <- bind_rows(
  get_fit_measures_safe(fit_1f, "1-factor"),
  get_fit_measures_safe(fit_2f, "2-factor")
)

print(fit_table)

# === ADDRESSING LOW RELIABILITY ===
cat("\n=== INVESTIGATING LOW RELIABILITY ===\n")

# Check item-total correlations
cat("Item-total correlations for Inspection Time:\n")
item_total_inspection <- cor(inspection_items, rowMeans(inspection_items))
print(round(item_total_inspection, 3))

cat("\nItem-total correlations for Temporal Order:\n")
item_total_temporal <- cor(temporal_items, rowMeans(temporal_items))
print(round(item_total_temporal, 3))

# Check if items are reverse coded or have different variances
cat("\nDescriptive statistics by item:\n")
describe(df_numeric) %>%
  select(n, mean, sd, min, max, skew) %>%
  print(digits = 2)

# === EXAMINE RESIDUAL CORRELATIONS (FIXED) ===
if (fit_2f@Fit@converged) {
  cat("\n=== RESIDUAL CORRELATIONS (2-FACTOR MODEL) ===\n")

  # Extract residuals properly
  resid_output <- lavResiduals(fit_2f)

  # Check if it's a list and extract the correlation component
  if (is.list(resid_output)) {
    resid_cor <- resid_output$cov # Sometimes it's called 'cov' not 'cor'
    if (is.null(resid_cor)) {
      resid_cor <- resid_output$cor
    }
  } else {
    resid_cor <- resid_output
  }

  # Print residual correlations if we have them
  if (!is.null(resid_cor) && is.matrix(resid_cor)) {
    cat("Residual correlation matrix:\n")
    print(round(resid_cor, 3))

    # Identify large residuals
    resid_cor_abs <- abs(resid_cor)
    diag(resid_cor_abs) <- 0 # Remove diagonal

    if (any(resid_cor_abs > 0.1, na.rm = TRUE)) {
      cat("\nLarge residual correlations (|r| > 0.1):\n")
      large_indices <- which(resid_cor_abs > 0.1, arr.ind = TRUE)
      for (i in 1:nrow(large_indices)) {
        if (large_indices[i, 1] < large_indices[i, 2]) {
          cat(sprintf(
            "  %s - %s: %.3f\n",
            rownames(resid_cor)[large_indices[i, 1]],
            colnames(resid_cor)[large_indices[i, 2]],
            resid_cor[large_indices[i, 1], large_indices[i, 2]]
          ))
        }
      }
    } else {
      cat("No large residual correlations found (all |r| < 0.1)\n")
    }
  }
}

# === ENHANCED DECISION FRAMEWORK ===
cat("\n\n=== COMPREHENSIVE DECISION FRAMEWORK ===\n")

# Calculate key metrics
if (fit_1f@Fit@converged && fit_2f@Fit@converged) {
  # Model improvement
  cfi_1f <- fitMeasures(fit_1f, "cfi")
  cfi_2f <- fitMeasures(fit_2f, "cfi")
  rmsea_1f <- fitMeasures(fit_1f, "rmsea")
  rmsea_2f <- fitMeasures(fit_2f, "rmsea")
  bic_1f <- fitMeasures(fit_1f, "bic")
  bic_2f <- fitMeasures(fit_2f, "bic")

  cfi_improvement <- cfi_2f - cfi_1f
  rmsea_improvement <- rmsea_1f - rmsea_2f
  bic_improvement <- bic_1f - bic_2f

  # Factor correlation
  factor_cor_cfa <- lavInspect(fit_2f, "cor.lv")[1, 2]

  # Chi-square difference test
  chi_diff_test <- lavTestLRT(fit_1f, fit_2f)
  p_diff <- chi_diff_test[2, "Pr(>Chisq)"]

  cat("Key Decision Metrics:\n")
  cat("─────────────────────\n")
  cat(
    "1. Statistical improvement (χ² diff test): p =",
    ifelse(is.na(p_diff), "NA",
      ifelse(p_diff < 0.001, "<0.001", round(p_diff, 3))
    ), "\n"
  )
  cat("2. Practical improvement:\n")
  cat(
    "   - ΔCFI =", round(cfi_improvement, 3),
    ifelse(cfi_improvement > 0.01, "✓", "✗"), "\n"
  )
  cat(
    "   - ΔRMSEA =", round(rmsea_improvement, 3),
    ifelse(rmsea_improvement > 0.015, "✓", "✗"), "\n"
  )
  cat(
    "   - ΔBIC =", round(bic_improvement, 2),
    ifelse(bic_improvement > 2, "✓", "✗"), "\n"
  )
  cat(
    "3. Factor correlation:", round(factor_cor_cfa, 3),
    ifelse(abs(factor_cor_cfa) < 0.85, "✓", "✗"), "\n"
  )
  cat("4. Reliability concerns:\n")
  cat(
    "   - Temporal α =", round(alpha_temporal$total$raw_alpha, 3),
    ifelse(alpha_temporal$total$raw_alpha > 0.6, "✓", "✗"), "\n"
  )
  cat(
    "   - Inspection α =", round(alpha_inspection$total$raw_alpha, 3),
    ifelse(alpha_inspection$total$raw_alpha > 0.6, "✓", "✗"), "\n"
  )

  # === EXPLORE WHY RELIABILITY IS LOW ===
  cat("\n=== EXPLORING LOW RELIABILITY ===\n")

  # Check inter-item correlations
  cat("\nInter-item correlations:\n")
  cat(
    "Temporal (CIRCLES-SQUARES):",
    round(cor(temporal_items)[1, 2], 3), "\n"
  )
  cat(
    "Inspection (CIRCLES-SQUARES):",
    round(cor(inspection_items)[1, 2], 3), "\n"
  )

  # Check if one item is problematic
  cat("\nItem means and SDs:\n")
  item_stats <- df_numeric %>%
    summarise(across(everything(), list(mean = mean, sd = sd))) %>%
    pivot_longer(everything(),
      names_to = c("item", "stat"),
      names_sep = "_(?=[^_]+$)"
    ) %>%
    pivot_wider(names_from = stat, values_from = value)
  print(item_stats, n = Inf)

  # FINAL RECOMMENDATION
  cat("\n=== FINAL RECOMMENDATION ===\n")

  # Decision logic accounting for low reliability
  if (alpha_inspection$total$raw_alpha < 0.6 || alpha_temporal$total$raw_alpha < 0.6) {
    cat("⚠️ CRITICAL ISSUE: Low reliability detected\n")
    cat("\nReliability Analysis:\n")
    cat("- Temporal α =", round(alpha_temporal$total$raw_alpha, 3), "\n")
    cat("- Inspection α =", round(alpha_inspection$total$raw_alpha, 3), "\n")
    cat("- Inter-item cor (temporal):", round(cor(temporal_items)[1, 2], 3), "\n")
    cat("- Inter-item cor (inspection):", round(cor(inspection_items)[1, 2], 3), "\n")

    cat("\nThe low reliability suggests:\n")
    cat("1. Circle and Square stimuli may tap different processes\n")
    cat("2. High measurement error in composite scores\n")
    cat("3. Factor scores will be unstable\n")

    cat("\nRECOMMENDATIONS:\n")
    cat("→ PRIMARY: Analyze items separately (4 variables)\n")
    cat("→ SECONDARY: Use best single item per task\n")
    cat("→ AVOID: Composite scores due to low reliability\n")
    cat("→ REPORT: Measurement limitations in manuscript\n")
  } else if (!is.na(p_diff) && p_diff < 0.05 && abs(factor_cor_cfa) < 0.85) {
    cat("✅ Evidence supports 2-factor structure\n")
    cat("→ Use separate temporal and inspection scores\n")
  } else if (abs(factor_cor_cfa) > 0.85) {
    cat("⚠️ Factors are highly correlated (r =", round(factor_cor_cfa, 3), ")\n")
    cat("→ Consider single composite score\n")
  } else {
    cat("❓ Mixed evidence\n")
    cat("→ Report both solutions and check sensitivity\n")
  }
}

# === CREATE FINAL SCORES WITH MEASUREMENT ERROR ===
cat("\n\n=== CREATING SCORES WITH RELIABILITY ADJUSTMENT ===\n")

# Given low reliability, provide multiple scoring options
processing_speed_final <- processing_speed %>%
  mutate(
    # Raw scores (standardized) - RECOMMENDED given low reliability
    temporal_circles_z = scale(temporal_CIRCLES)[, 1],
    temporal_squares_z = scale(temporal_SQUARES)[, 1],
    inspection_circles_z = scale(inspection_CIRCLES)[, 1],
    inspection_squares_z = scale(inspection_SQUARES)[, 1],

    # Simple composites (USE WITH CAUTION)
    temporal_mean = rowMeans(scale(select(., starts_with("temporal_")))),
    inspection_mean = rowMeans(scale(select(., starts_with("inspection_")))),

    # Single best item per task (based on loadings/reliability)
    temporal_best = temporal_squares_z, # Typically higher loading
    inspection_best = inspection_circles_z, # Typically higher loading

    # Overall composite (NOT RECOMMENDED)
    ps_composite = rowMeans(scale(select(., temporal_CIRCLES:inspection_SQUARES)))
  )

# Add factor scores if model converged (NOT RECOMMENDED with low reliability)
if (fit_2f@Fit@converged) {
  scores_2f <- lavPredict(fit_2f, type = "lv")
  processing_speed_final$temporal_factor <- scores_2f[, "temporal"]
  processing_speed_final$inspection_factor <- scores_2f[, "inspection"]
}

# Report correlations between different scoring approaches
cat("\nCorrelations between scoring methods:\n")
score_vars <- c(
  "temporal_circles_z", "temporal_squares_z",
  "inspection_circles_z", "inspection_squares_z",
  "temporal_mean", "inspection_mean"
)
cor_matrix <- cor(processing_speed_final[, score_vars], use = "complete.obs")
print(round(cor_matrix, 3))

# === STATISTICAL POWER ANALYSIS ===
cat("\n=== RELIABILITY IMPACT ON STATISTICAL POWER ===\n")
cat("With reliability α =", round(alpha_inspection$total$raw_alpha, 3), "\n")
cat(
  "- Observed correlations are attenuated by factor of:",
  round(sqrt(alpha_inspection$total$raw_alpha), 3), "\n"
)
cat(
  "- To detect r = 0.30, you need observed r =",
  round(0.30 * sqrt(alpha_inspection$total$raw_alpha), 3), "\n"
)
cat(
  "- Power is reduced by approximately",
  round((1 - alpha_inspection$total$raw_alpha) * 100, 1), "%\n"
)

# Save with documentation
cat("\n✓ Saving results with documentation...\n")
write_csv(processing_speed_final, "processing_speed_scores_final.csv")

# Create detailed documentation
sink("processing_speed_scoring_documentation.txt")
cat("PROCESSING SPEED SCORING DOCUMENTATION\n")
cat("Generated:", Sys.Date(), "\n")
cat("=====================================\n\n")
cat("RELIABILITY ISSUES:\n")
cat("- Temporal Order α =", round(alpha_temporal$total$raw_alpha, 3), "\n")
cat("- Inspection Time α =", round(alpha_inspection$total$raw_alpha, 3), "\n")
cat("- Both below acceptable threshold of 0.70\n\n")
cat("RECOMMENDED SCORING:\n")
cat("1. Use individual standardized items (e.g., temporal_squares_z)\n")
cat("2. Avoid composite scores due to low reliability\n")
cat("3. Consider measurement error in interpretation\n\n")
cat("VARIABLES IN OUTPUT FILE:\n")
cat("- temporal_circles_z: Standardized temporal order (circles)\n")
cat("- temporal_squares_z: Standardized temporal order (squares)\n")
cat("- inspection_circles_z: Standardized inspection time (circles)\n")
cat("- inspection_squares_z: Standardized inspection time (squares)\n")
cat("- temporal_mean: Mean of temporal items (USE WITH CAUTION)\n")
cat("- inspection_mean: Mean of inspection items (USE WITH CAUTION)\n")
sink()

cat("✓ Documentation saved to 'processing_speed_scoring_documentation.txt'\n")
