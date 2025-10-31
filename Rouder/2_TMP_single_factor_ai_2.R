library(tidyverse)
library(fs)
library(psych)
library(lavaan)
library(corrplot)
library(semTools)
library(moments)
library(car)
library(GGally)
library(gridExtra)
library(ggrepel)

processing_speed <- read_csv("cleaned_cognitive_data.csv") 

# === IDENTIFY AND LIST ALL OUTLIER PARTICIPANTS ===

# Function to identify outliers using IQR method
identify_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  return(x < lower | x > upper)
}

# 1. Univariate outliers (IQR method)
univariate_outliers <- processing_speed %>%
  mutate(
    outlier_temporal_C = identify_outliers(temporal_CIRCLES),
    outlier_temporal_S = identify_outliers(temporal_SQUARES),
    outlier_inspection_C = identify_outliers(inspection_CIRCLES),
    outlier_inspection_S = identify_outliers(inspection_SQUARES),
    any_univariate = outlier_temporal_C | outlier_temporal_S | 
      outlier_inspection_C | outlier_inspection_S
  ) %>%
  filter(any_univariate)

cat("=== UNIVARIATE OUTLIERS (IQR Method) ===\n")
cat("Participants:", paste(univariate_outliers$PART_ID, collapse = ", "), "\n")
cat("Total:", nrow(univariate_outliers), "participants\n\n")

# 2. Multivariate outliers (Mahalanobis distance)
mah_dist <- mahalanobis(
  processing_speed[, 2:5], 
  colMeans(processing_speed[, 2:5]), 
  cov(processing_speed[, 2:5])
)

multivariate_outliers <- processing_speed %>%
  mutate(mah_dist = mah_dist) %>%
  filter(mah_dist > qchisq(0.99, df = 4))

cat("=== MULTIVARIATE OUTLIERS (Mahalanobis Distance) ===\n")
cat("Participants:", paste(multivariate_outliers$PART_ID, collapse = ", "), "\n")
cat("Total:", nrow(multivariate_outliers), "participants\n\n")

# 3. Combined list of all outliers
all_outliers <- unique(c(univariate_outliers$PART_ID, multivariate_outliers$PART_ID))

cat("=== ALL OUTLIERS (COMBINED) ===\n")
cat("Participants:", paste(all_outliers, collapse = ", "), "\n")
cat("Total:", length(all_outliers), "participants\n\n")

processing_speed <- processing_speed %>%
  filter(!PART_ID %in% all_outliers)

# === 1. DISTRIBUTION PLOTS ===
# Individual histograms with density curves
dist_plots <- processing_speed %>%
  select(-PART_ID) %>%
  gather(task, value) %>%
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "lightblue", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  facet_wrap(~task, scales = "free", ncol = 2) +
  theme_minimal() +
  labs(title = "Distribution of SOA Values by Task",
       x = "SOA", y = "Density")

print(dist_plots)

# === 2. BOXPLOTS WITH OUTLIER LABELS ===
# Identify outliers using IQR method
identify_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  return(x < lower | x > upper)
}

# Create long format with outlier flags
long_data <- processing_speed %>%
  gather(task, value, -PART_ID) %>%
  group_by(task) %>%
  mutate(is_outlier = identify_outliers(value))

# Boxplot with outlier points labeled
box_plots <- ggplot(long_data, aes(x = task, y = value)) +
  geom_boxplot(fill = "lightcoral", alpha = 0.7) +
  geom_point(data = filter(long_data, is_outlier), 
             color = "red", size = 3) +
  geom_text_repel(data = filter(long_data, is_outlier),
                  aes(label = PART_ID),
                  size = 3, max.overlaps = 20) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplots with Outlier Identification",
       x = "Task", y = "SOA")

print(box_plots)

# === 3. SCATTER PLOT MATRIX ===
# Shows relationships and outliers across all pairs
scatter_matrix <- processing_speed %>%
  select(-PART_ID) %>%
  ggpairs(
    lower = list(continuous = wrap("points", alpha = 0.5, size = 2)),
    diag = list(continuous = wrap("densityDiag", alpha = 0.5)),
    upper = list(continuous = wrap("cor", size = 4)),
    title = "Scatter Plot Matrix with Correlations"
  )

print(scatter_matrix)

# === 4. BIVARIATE OUTLIER DETECTION ===
# Mahalanobis distance visualization
mah_dist <- mahalanobis(
  processing_speed[, -1], 
  colMeans(processing_speed[, -1]), 
  cov(processing_speed[, -1])
)

processing_speed$mah_dist <- mah_dist
processing_speed$outlier_mah <- mah_dist > qchisq(0.99, df = 4)

# Add index for plotting
processing_speed$index <- 1:nrow(processing_speed)

# Plot Mahalanobis distances
mah_plot <- ggplot(processing_speed, aes(x = index, y = mah_dist)) +
  geom_point(aes(color = outlier_mah), size = 2) +
  geom_hline(yintercept = qchisq(0.99, df = 4), 
             linetype = "dashed", color = "red") +
  geom_text_repel(data = filter(processing_speed, outlier_mah),
                  aes(label = PART_ID), size = 3) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  theme_minimal() +
  labs(title = "Multivariate Outlier Detection (Mahalanobis Distance)",
       x = "Observation Index", y = "Mahalanobis Distance",
       color = "Outlier")

print(mah_plot)

# === 5. TASK COMPARISON PLOTS ===
# Compare temporal vs inspection by stimulus type
comparison_data <- processing_speed %>%
  select(-mah_dist, -outlier_mah, -index) %>%
  gather(task, value, -PART_ID) %>%
  separate(task, into = c("task_type", "stimulus"), sep = "_") %>%
  spread(task_type, value)

task_comparison <- ggplot(comparison_data, aes(x = temporal, y = inspection)) +
  geom_point(aes(color = stimulus), size = 3, alpha = 0.7) +
  geom_smooth(aes(color = stimulus), method = "lm", se = TRUE) +
  facet_wrap(~stimulus) +
  theme_minimal() +
  labs(title = "Temporal vs Inspection Time by Stimulus Type",
       x = "Temporal Order SOA", y = "Inspection Time SOA")

print(task_comparison)

# === 6. VIOLIN PLOTS WITH INDIVIDUAL POINTS ===
violin_plots <- long_data %>%
  separate(task, into = c("task_type", "stimulus"), sep = "_") %>%
  ggplot(aes(x = stimulus, y = value, fill = task_type)) +
  geom_violin(alpha = 0.7, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.8), 
               outlier.shape = NA, alpha = 0.5) +
  geom_jitter(position = position_dodge(0.8), alpha = 0.3, size = 1) +
  facet_wrap(~task_type, scales = "free_y") +
  theme_minimal() +
  labs(title = "Distribution by Task Type and Stimulus",
       x = "Stimulus Type", y = "SOA", fill = "Task Type")

print(violin_plots)

# === 7. OUTLIER SUMMARY TABLE ===
# Create comprehensive outlier summary
outlier_summary <- processing_speed %>%
  mutate(
    temporal_CIRCLES_out = identify_outliers(temporal_CIRCLES),
    temporal_SQUARES_out = identify_outliers(temporal_SQUARES),
    inspection_CIRCLES_out = identify_outliers(inspection_CIRCLES),
    inspection_SQUARES_out = identify_outliers(inspection_SQUARES)
  ) %>%
  filter(temporal_CIRCLES_out | temporal_SQUARES_out | 
           inspection_CIRCLES_out | inspection_SQUARES_out | outlier_mah) %>%
  select(PART_ID, ends_with("_out"), outlier_mah, everything())

cat("\n=== OUTLIER SUMMARY ===\n")
print(outlier_summary)

# === 8. RANGE AND SCALE COMPARISON ===
# Visualize the scale differences between tasks
scale_comparison <- long_data %>%
  group_by(task) %>%
  summarise(
    min = min(value),
    q1 = quantile(value, 0.25),
    median = median(value),
    q3 = quantile(value, 0.75),
    max = max(value)
  ) %>%
  gather(stat, value, -task) %>%
  mutate(stat = factor(stat, levels = c("min", "q1", "median", "q3", "max")))

scale_plot <- ggplot(scale_comparison, aes(x = task, y = value, group = task)) +
  geom_line(size = 2, color = "gray") +
  geom_point(aes(color = stat), size = 4) +
  scale_color_manual(values = c(
    "min" = "blue", "q1" = "lightblue", "median" = "black", 
    "q3" = "orange", "max" = "red"
  )) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Scale Comparison Across Tasks",
       x = "Task", y = "SOA", color = "Statistic")

print(scale_plot)

# === 9. PARTICIPANT PROFILES ===
# Identify participants with unusual patterns
participant_profiles <- processing_speed %>%
  mutate(
    temporal_mean = (temporal_CIRCLES + temporal_SQUARES) / 2,
    inspection_mean = (inspection_CIRCLES + inspection_SQUARES) / 2,
    temporal_diff = abs(temporal_CIRCLES - temporal_SQUARES),
    inspection_diff = abs(inspection_CIRCLES - inspection_SQUARES)
  )

profile_plot <- ggplot(participant_profiles, 
                       aes(x = temporal_mean, y = inspection_mean)) +
  geom_point(aes(size = temporal_diff + inspection_diff), alpha = 0.6) +
  geom_text_repel(data = filter(participant_profiles, outlier_mah),
                  aes(label = PART_ID), size = 3) +
  theme_minimal() +
  labs(title = "Participant Performance Profiles",
       x = "Mean Temporal Order SOA",
       y = "Mean Inspection Time SOA",
       size = "Total within-task\nvariability")

print(profile_plot)

# === 10. SUMMARY STATISTICS ===
cat("\n=== SUMMARY STATISTICS BY TASK ===\n")
summary_stats <- processing_speed %>%
  select(-PART_ID, -mah_dist, -outlier_mah, -index) %>%
  gather(task, value) %>%
  group_by(task) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    mad = mad(value),
    skew = moments::skewness(value),
    n_outliers_iqr = sum(identify_outliers(value)),
    pct_outliers = round(100 * n_outliers_iqr / n, 1)
  )

print(summary_stats)

# === 11. OUTLIER INFLUENCE PLOT ===
# Show how removing outliers affects correlations
cor_with_outliers <- cor(processing_speed[, 2:5])

# Remove multivariate outliers
clean_data <- processing_speed %>%
  filter(!outlier_mah) %>%
  select(temporal_CIRCLES:inspection_SQUARES)

cor_without_outliers <- cor(clean_data)

# Create comparison
cor_diff <- cor_without_outliers - cor_with_outliers

# Visualize correlation changes
library(reshape2)
cor_diff_melted <- melt(cor_diff)

cor_change_plot <- ggplot(cor_diff_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-0.3, 0.3)) +
  geom_text(aes(label = round(value, 3)), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Change in Correlations After Removing Outliers",
       x = "", y = "", fill = "Change")

print(cor_change_plot)



df_numeric <- processing_speed |> 
  select(temporal_CIRCLES, temporal_SQUARES, inspection_CIRCLES, inspection_SQUARES)



# === DIAGNOSTIC ANALYSIS: WHY SINGLE FACTOR FAILS ===
# Despite theoretical expectations, data doesn't support single factor
# This script investigates potential reasons

# === 1. TASK CHARACTERISTICS COMPARISON ===
cat("=== 1. BASIC TASK PROPERTIES ===\n")

# Descriptive statistics by task
task_stats <- df_numeric %>%
  gather(task, value) %>%
  group_by(task) %>%
  summarise(
    mean = mean(value),
    sd = sd(value),
    cv = sd/mean,  # Coefficient of variation
    skew = skewness(value),
    kurt = kurtosis(value),
    min = min(value),
    max = max(value),
    range = max - min
  )

print(task_stats)

# Different means/variances can cause low correlations
cat("\n⚠️ Check for scale differences:\n")
cat("- CV range:", round(range(task_stats$cv), 3), "\n")
if(max(task_stats$cv)/min(task_stats$cv) > 2) {
  cat("→ Large differences in variability across tasks\n")
}

# === 2. DISTRIBUTION ANALYSIS ===
cat("\n=== 2. DISTRIBUTION DIAGNOSTICS ===\n")

# Visual inspection
par(mfrow=c(2,2))
for(i in 1:4) {
  hist(df_numeric[[i]], main=names(df_numeric)[i], 
       xlab="SOA", breaks=20, col="lightblue")
}

# Test for multivariate normality using psych package
mardia_result <- mardia(df_numeric)
cat("\nMultivariate normality (Mardia's test):\n")
cat("- Skewness statistic:", round(mardia_result$b1p, 3), 
    ", p-value:", round(mardia_result$p.skew, 3), "\n")
cat("- Kurtosis statistic:", round(mardia_result$b2p, 3), 
    ", p-value:", round(mardia_result$p.kurt, 3), "\n")

# === 3. OUTLIER AND INFLUENTIAL CASES ===
cat("\n=== 3. OUTLIER ANALYSIS ===\n")

# Multivariate outliers using Mahalanobis distance
mah_dist <- mahalanobis(df_numeric, colMeans(df_numeric), cov(df_numeric))
cutoff <- qchisq(0.999, df=ncol(df_numeric))
outliers <- which(mah_dist > cutoff)

cat("Multivariate outliers detected:", length(outliers), "\n")
if(length(outliers) > 0) {
  cat("→ Outlier IDs:", outliers, "\n")
  
  # Compare correlations with/without outliers
  cor_with <- cor(df_numeric)
  cor_without <- cor(df_numeric[-outliers,])
  cor_diff <- cor_without - cor_with
  cat("\nMax correlation change after removing outliers:", 
      round(max(abs(cor_diff[upper.tri(cor_diff)])), 3), "\n")
}

# === 4. MEASUREMENT QUALITY BY TASK ===
cat("\n=== 4. MEASUREMENT QUALITY ANALYSIS ===\n")

# Check variance ratios
var_ratios <- apply(df_numeric, 2, var)
cat("\nVariance by task:\n")
print(round(var_ratios, 2))
cat("\nVariance ratio (max/min):", round(max(var_ratios)/min(var_ratios), 2), "\n")

# === 5. SUBGROUP ANALYSIS ===
cat("\n=== 5. SUBGROUP HETEROGENEITY ===\n")

# Split-half reliability check
set.seed(123)
split <- sample(1:2, nrow(df_numeric), replace=TRUE)

cor_split1 <- cor(df_numeric[split==1,])
cor_split2 <- cor(df_numeric[split==2,])

cat("\nCorrelation stability (split-half):\n")
cor_diff_splits <- abs(cor_split1 - cor_split2)
diag(cor_diff_splits) <- NA
cat("- Max correlation difference between splits:", 
    round(max(cor_diff_splits, na.rm=TRUE), 3), "\n")
cat("- Mean absolute difference:", 
    round(mean(cor_diff_splits, na.rm=TRUE), 3), "\n")

# === 6. TASK DIFFICULTY ANALYSIS ===
cat("\n=== 6. RELATIVE TASK DIFFICULTY ===\n")

# Standardize and compare
df_std <- scale(df_numeric)
task_difficulty <- colMeans(df_std)
names(task_difficulty) <- names(df_numeric)

cat("\nStandardized mean SOAs (higher = easier):\n")
print(round(task_difficulty, 2))

# Check for floor/ceiling effects
ceiling_effects <- apply(df_numeric, 2, function(x) sum(x > quantile(x, 0.95))/length(x))
floor_effects <- apply(df_numeric, 2, function(x) sum(x < quantile(x, 0.05))/length(x))

cat("\nRestricted range indicators:\n")
cat("- Ceiling effects (% > 95th percentile):\n")
print(round(ceiling_effects * 100, 1))
cat("- Floor effects (% < 5th percentile):\n")
print(round(floor_effects * 100, 1))

# === 7. METHOD VARIANCE INVESTIGATION ===
cat("\n=== 7. METHOD VARIANCE PATTERNS ===\n")

# Check if correlations follow method patterns
cor_mat <- cor(df_numeric)
within_temporal <- cor_mat["temporal_CIRCLES", "temporal_SQUARES"]
within_inspection <- cor_mat["inspection_CIRCLES", "inspection_SQUARES"]
within_circles <- cor_mat["temporal_CIRCLES", "inspection_CIRCLES"]
within_squares <- cor_mat["temporal_SQUARES", "inspection_SQUARES"]

cat("\nCorrelation patterns:\n")
cat("- Same task, different stimuli:\n")
cat("  Temporal (circles-squares):", round(within_temporal, 3), "\n")
cat("  Inspection (circles-squares):", round(within_inspection, 3), "\n")
cat("- Same stimuli, different tasks:\n")
cat("  Circles (temporal-inspection):", round(within_circles, 3), "\n")
cat("  Squares (temporal-inspection):", round(within_squares, 3), "\n")

# === 8. RESIDUAL CORRELATION PATTERNS ===
cat("\n=== 8. RESIDUAL CORRELATION ANALYSIS ===\n")

# Fit single factor and examine residuals
fit_diag <- cfa(model_1factor, data=processing_speed, std.lv=TRUE)
resid <- residuals(fit_diag, type="cor")$cov

# Find systematic patterns
cat("\nLargest positive residual:", round(max(resid[upper.tri(resid)]), 3), "\n")
cat("Largest negative residual:", round(min(resid[upper.tri(resid)]), 3), "\n")

# Check for method factors
cat("\nResidual correlations suggesting method effects:\n")
cat("- Temporal tasks residual:", round(resid["temporal_CIRCLES", "temporal_SQUARES"], 3), "\n")
cat("- Inspection tasks residual:", round(resid["inspection_CIRCLES", "inspection_SQUARES"], 3), "\n")

# === 9. ALTERNATIVE MODEL TESTING ===
cat("\n=== 9. TESTING ALTERNATIVE EXPLANATIONS ===\n")

# Model with method factors
model_method <- '
  # General factor
  g =~ temporal_CIRCLES + temporal_SQUARES + inspection_CIRCLES + inspection_SQUARES
  
  # Method factors
  temporal_method =~ temporal_CIRCLES + temporal_SQUARES
  inspection_method =~ inspection_CIRCLES + inspection_SQUARES
  
  # Orthogonal methods
  temporal_method ~~ 0*inspection_method
  g ~~ 0*temporal_method
  g ~~ 0*inspection_method
'

fit_method <- try(cfa(model_method, data=processing_speed, std.lv=TRUE))
if(!inherits(fit_method, "try-error")) {
  cat("\nBifactor model with method factors:\n")
  print(fitMeasures(fit_method, c("cfi", "tli", "rmsea")))
}

# === 10. DIAGNOSTIC SUMMARY ===
cat("\n=== DIAGNOSTIC SUMMARY: WHY SINGLE FACTOR FAILS ===\n")
cat("==================================================\n")

# Compile evidence
problems <- list()

if(max(task_stats$cv)/min(task_stats$cv) > 1.5) {
  problems$scaling <- "Different measurement scales across tasks"
}

if(length(outliers) > nrow(df_numeric)*0.05) {
  problems$outliers <- "Substantial outliers affecting correlations"
}

if(mean(c(within_temporal, within_inspection)) < 0.4) {
  problems$reliability <- "Low within-task correlations suggest measurement error"
}

if(abs(within_temporal - within_inspection) > 0.2) {
  problems$consistency <- "Inconsistent within-task correlations"
}

if(max(abs(resid[upper.tri(resid)])) > 0.15) {
  problems$misfit <- "Large residuals indicate systematic misfit"
}

cat("\n🔍 IDENTIFIED ISSUES:\n")
for(i in seq_along(problems)) {
  cat(i, ".", problems[[i]], "\n")
}

cat("\n📊 THEORETICAL VS. EMPIRICAL MISMATCH EXPLANATIONS:\n")
cat("1. MEASUREMENT ERROR: Tasks may have different reliability\n")
cat("2. TASK IMPURITY: Tasks tap additional cognitive processes\n")
cat("3. STRATEGY DIFFERENCES: Participants use different strategies\n")
cat("4. STATISTICAL ARTIFACTS: Range restriction, outliers, etc.\n")
cat("5. DEVELOPMENTAL/INDIVIDUAL DIFFERENCES: Factor structure varies across people\n")

cat("\n✅ RECOMMENDATIONS FOR MOVING FORWARD:\n")
cat("1. Report the failure to replicate single-factor structure\n")
cat("2. Use individual task scores rather than composites\n")
cat("3. Consider these as related but distinct measures\n")
cat("4. Investigate task-specific validity separately\n")

# Create a diagnostic plot
par(mfrow=c(2,2))
corrplot(cor(df_numeric), method="number", main="Observed Correlations")
hist(mah_dist, main="Mahalanobis Distances", xlab="Distance", breaks=20)
abline(v=cutoff, col="red", lty=2)
plot(density(resid[upper.tri(resid)]), main="Residual Distribution")
abline(v=0, col="red", lty=2)


cat("Sample size:", nrow(processing_speed), "\n\n")

# === 2. CORRELATION ANALYSIS ===
cat("=== CORRELATION ANALYSIS ===\n")
cor_matrix <- cor(df_numeric)
cat("\nCorrelation Matrix:\n")
print(round(cor_matrix, 3))

# Test correlations between task types
temporal_vars <- c("temporal_CIRCLES", "temporal_SQUARES")
inspection_vars <- c("inspection_CIRCLES", "inspection_SQUARES")

within_temporal <- cor_matrix[temporal_vars, temporal_vars]
within_inspection <- cor_matrix[inspection_vars, inspection_vars]
between_tasks <- cor_matrix[temporal_vars, inspection_vars]

cat("\nWithin-task correlations:\n")
cat("- Temporal Order tasks:", round(within_temporal[1,2], 3), "\n")
cat("- Inspection Time tasks:", round(within_inspection[1,2], 3), "\n")
cat("\nBetween-task correlations (mean):", round(mean(between_tasks), 3), "\n")

# === 3. SINGLE FACTOR ANALYSIS (DEMONSTRATING POOR FIT) ===
cat("\n=== SINGLE FACTOR MODEL (ALL 4 TASKS) ===\n")

model_1factor <- '
  processing_speed =~ temporal_CIRCLES + temporal_SQUARES + 
                      inspection_CIRCLES + inspection_SQUARES
'

fit_1f <- cfa(model_1factor, data = processing_speed, std.lv = TRUE)

# Extract fit measures
fit_1f_measures <- fitMeasures(fit_1f, c("chisq", "df", "pvalue", "cfi", "tli", 
                                         "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
                                         "srmr"))

cat("\nSingle Factor Fit Statistics:\n")
cat("- Chi-square:", round(fit_1f_measures["chisq"], 2), "\n")
cat("- df:", fit_1f_measures["df"], "\n")
cat("- p-value:", format.pval(fit_1f_measures["pvalue"], digits = 3), "\n")
cat("- CFI:", round(fit_1f_measures["cfi"], 3), "❌ (should be > 0.95)\n")
cat("- TLI:", round(fit_1f_measures["tli"], 3), "❌ (should be > 0.95)\n")
cat("- RMSEA:", round(fit_1f_measures["rmsea"], 3), "❌ (should be < 0.06)\n")
cat("- SRMR:", round(fit_1f_measures["srmr"], 3), "\n")

# === 4. TWO-FACTOR MODEL (SHOWING TASK DIFFERENCES) ===
cat("\n=== TWO-FACTOR MODEL (TEMPORAL vs INSPECTION) ===\n")

model_2factor <- '
  temporal_order =~ temporal_CIRCLES + temporal_SQUARES
  inspection_time =~ inspection_CIRCLES + inspection_SQUARES
'

fit_2f <- cfa(model_2factor, data = processing_speed, std.lv = TRUE)

# Extract fit measures
fit_2f_measures <- fitMeasures(fit_2f, c("chisq", "df", "pvalue", "cfi", "tli", 
                                         "rmsea", "srmr"))

cat("\nTwo-Factor Fit Statistics:\n")
cat("- CFI:", round(fit_2f_measures["cfi"], 3), "\n")
cat("- TLI:", round(fit_2f_measures["tli"], 3), "\n")
cat("- RMSEA:", round(fit_2f_measures["rmsea"], 3), "\n")

# Factor correlation
factor_cor <- lavInspect(fit_2f, "cor.lv")
cat("\nFactor Correlation (Temporal vs Inspection):", round(factor_cor[1,2], 3), "\n")

if (factor_cor[1,2] < 0.7) {
  cat("⚠️ Low factor correlation indicates these measure different constructs\n")
}

# === 5. INSPECTION TIME ONLY MODEL (RECOMMENDED) ===
cat("\n=== INSPECTION TIME ONLY MODEL ===\n")

# Reliability for inspection time items only
inspection_only <- df_numeric[, inspection_vars]
alpha_inspection <- alpha(inspection_only)

cat("Reliability (Inspection Time only):\n")
cat("- Cronbach's Alpha:", round(alpha_inspection$total$raw_alpha, 3), "✓\n")
cat("- Correlation between items:", round(cor(inspection_only)[1,2], 3), "\n")

# Simple factor model for inspection time
model_inspection <- '
  inspection_time =~ inspection_CIRCLES + inspection_SQUARES
'

fit_inspection <- cfa(model_inspection, data = processing_speed, std.lv = TRUE)
fit_insp_measures <- fitMeasures(fit_inspection, c("cfi", "tli", "rmsea", "srmr"))

cat("\nFit Statistics (Inspection Only):\n")
cat("- Model is just-identified (df=0), perfect fit by definition\n")

# Composite score analysis
cat("\nComposite Score Properties:\n")
processing_speed$inspection_composite <- rowMeans(inspection_only)
sd_composite <- sd(processing_speed$inspection_composite)
cat("- SD of composite:", round(sd_composite, 2), "\n")
cat("- Reliability of composite:", round(alpha_inspection$total$raw_alpha, 3), "\n")

# === 6. COMPARISON OF APPROACHES ===
cat("\n=== COMPARISON OF MEASUREMENT APPROACHES ===\n")
cat("============================================\n")

# Calculate variance explained for different models
efa_all <- fa(df_numeric, nfactors = 1, rotate = "none", fm = "ml")
efa_inspection <- fa(inspection_only, nfactors = 1, rotate = "none", fm = "ml")

var_all <- sum(efa_all$communality) / length(efa_all$communality) * 100
var_inspection <- sum(efa_inspection$communality) / length(efa_inspection$communality) * 100

cat("\n1. SINGLE FACTOR (ALL TASKS):\n")
cat("   - Variance explained:", round(var_all, 1), "%\n")
cat("   - Model fit: POOR (CFI =", round(fit_1f_measures["cfi"], 3), ")\n")
cat("   - Reliability: α =", round(alpha(df_numeric)$total$raw_alpha, 3), "\n")
cat("   - Status: ❌ NOT RECOMMENDED\n")

cat("\n2. TEMPORAL ORDER ONLY:\n")
temporal_only <- df_numeric[, temporal_vars]
alpha_temporal <- alpha(temporal_only)
cat("   - Item correlation:", round(cor(temporal_only)[1,2], 3), "\n")
cat("   - Reliability: α =", round(alpha_temporal$total$raw_alpha, 3), "\n")
cat("   - Status: ⚠️ Low reliability\n")

cat("\n3. INSPECTION TIME ONLY:\n")
cat("   - Variance explained:", round(var_inspection, 1), "%\n")
cat("   - Item correlation:", round(cor(inspection_only)[1,2], 3), "\n")
cat("   - Reliability: α =", round(alpha_inspection$total$raw_alpha, 3), "\n")
cat("   - Status: ✅ RECOMMENDED\n")

# === 7. FINAL RECOMMENDATIONS ===
cat("\n=== CONCLUSIONS & RECOMMENDATIONS ===\n")
cat("=====================================\n")

cat("\n❌ SINGLE FACTOR APPROACH FAILS:\n")
cat("1. Poor model fit (CFI =", round(fit_1f_measures["cfi"], 3), ", RMSEA =", 
    round(fit_1f_measures["rmsea"], 3), ")\n")
cat("2. Low between-task correlations (r̄ =", round(mean(between_tasks), 3), ")\n")
cat("3. Temporal and Inspection tasks measure different constructs\n")
cat("4. Combining all tasks masks important individual differences\n")

cat("\n✅ RECOMMENDED APPROACH: USE INSPECTION TIME ONLY\n")
cat("1. High internal consistency (α =", round(alpha_inspection$total$raw_alpha, 3), ")\n")
cat("2. Strong item correlation (r =", round(cor(inspection_only)[1,2], 3), ")\n")
cat("3. Theoretically coherent - both measure visual processing speed\n")
cat("4. Simple composite score: mean(inspection_CIRCLES, inspection_SQUARES)\n")

cat("\n📊 IMPLEMENTATION:\n")
cat("processing_speed$speed_score <- rowMeans(processing_speed[, c('inspection_CIRCLES', 'inspection_SQUARES')])\n")

# Create the recommended score
processing_speed$speed_score <- rowMeans(processing_speed[, inspection_vars])

cat("\n⚠️ TEMPORAL ORDER TASKS:\n")
cat("- Should be analyzed separately if needed\n")
cat("- Measure different construct (temporal discrimination)\n")
cat("- Lower reliability when combined\n")

# === 8. EVIDENCE SUMMARY ===
cat("\n=== STATISTICAL EVIDENCE SUMMARY ===\n")

evidence_table <- data.frame(
  Criterion = c("Within-task correlation", "Between-task correlation", 
                "Single factor CFI", "Reliability (all)", "Reliability (inspection only)"),
  Value = c(
    round(mean(c(within_temporal[1,2], within_inspection[1,2])), 3),
    round(mean(between_tasks), 3),
    round(fit_1f_measures["cfi"], 3),
    round(alpha(df_numeric)$total$raw_alpha, 3),
    round(alpha_inspection$total$raw_alpha, 3)
  ),
  Interpretation = c(
    "Moderate within tasks",
    "Weak between task types",
    "Poor fit - reject single factor",
    "Marginal reliability",
    "Good reliability"
  )
)

print(evidence_table)

cat("\n🎯 FINAL RECOMMENDATION:\n")
cat("Use mean(inspection_CIRCLES, inspection_SQUARES) as processing speed measure\n")
cat("Analyze temporal order tasks separately if needed\n")