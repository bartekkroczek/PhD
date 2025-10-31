library(tidyverse)
library(broom)
library(car)
library(lmtest)
library(robustbase)
library(boot)
library(patchwork)

# Load data
rou_rep <- read_csv("saccades_factor.csv")

# 1. DESCRIPTIVE STATISTICS
cat("=== DESCRIPTIVE STATISTICS ===\n")
summary_stats <- rou_rep %>%
  summarise(
    across(c(inspection_time, anti), 
           list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE),
                min = ~min(.x, na.rm = TRUE),
                max = ~max(.x, na.rm = TRUE)), 
           .names = "{.col}_{.fn}")
  )
print(summary_stats)

# 2. ASSUMPTION CHECKS
cat("\n=== ASSUMPTION CHECKS ===\n")

# Normality (Shapiro-Wilk test)
shapiro_inspection <- shapiro.test(rou_rep$inspection_time)
shapiro_anti <- shapiro.test(rou_rep$anti)

cat("Normality Tests:\n")
cat("Inspection Time: W =", round(shapiro_inspection$statistic, 3), 
    ", p =", format.pval(shapiro_inspection$p.value, digits = 3), "\n")
cat("Anti: W =", round(shapiro_anti$statistic, 3), 
    ", p =", format.pval(shapiro_anti$p.value, digits = 3), "\n")

# Fit initial model for diagnostics
model_ols <- lm(anti ~ inspection_time, data = rou_rep)

# Heteroscedasticity (Breusch-Pagan test)
bp_test <- bptest(model_ols)
cat("\nHeteroscedasticity Test: BP =", round(bp_test$statistic, 3), 
    ", p =", format.pval(bp_test$p.value, digits = 3), "\n")

# Outliers (Bonferroni outlier test)
outliers <- outlierTest(model_ols)
cat("\nOutlier Test:\n")
print(outliers)

# 3. CORRELATIONS
cat("\n=== CORRELATIONS ===\n")
pearson_cor <- cor.test(rou_rep$inspection_time, rou_rep$anti, method = "pearson")
spearman_cor <- cor.test(rou_rep$inspection_time, rou_rep$anti, method = "spearman")

cat("Pearson: r =", round(pearson_cor$estimate, 3), 
    ", p =", format.pval(pearson_cor$p.value, digits = 3), "\n")
cat("Spearman: rho =", round(spearman_cor$estimate, 3), 
    ", p =", format.pval(spearman_cor$p.value, digits = 3), "\n")

# 4. REGRESSION MODELS
cat("\n=== REGRESSION MODELS ===\n")

# OLS regression
cat("OLS Regression:\n")
ols_summary <- summary(model_ols)
print(ols_summary)

# Robust regression
model_robust <- lmrob(anti ~ inspection_time, data = rou_rep)
cat("\nRobust Regression:\n")
robust_summary <- summary(model_robust)
print(robust_summary)

# 5. BOOTSTRAP CONFIDENCE INTERVALS
cat("\n=== BOOTSTRAP ANALYSIS ===\n")
set.seed(123)

# Bootstrap slope estimate
boot_slope <- function(data, indices) {
  coef(lm(anti ~ inspection_time, data = data[indices, ]))[2]
}

boot_results <- boot(rou_rep, boot_slope, R = 2000)
boot_ci <- boot.ci(boot_results, type = "bca")

cat("Bootstrap slope:", round(boot_results$t0, 3), "\n")
cat("Bootstrap 95% CI: [", round(boot_ci$bca[4], 3), ", ", 
    round(boot_ci$bca[5], 3), "]\n")

# Bootstrap p-value (permutation test)
boot_pvalue <- function(data, indices) {
  # Permute predictor to break association
  perm_data <- data
  perm_data$inspection_time <- sample(data$inspection_time)
  coef(lm(anti ~ inspection_time, data = perm_data))[2]
}

boot_null <- boot(rou_rep, boot_pvalue, R = 2000)
p_bootstrap <- mean(abs(boot_null$t) >= abs(boot_results$t0))
cat("Bootstrap p-value:", round(p_bootstrap, 3), "\n")

# 6. MODEL COMPARISON
cat("\n=== MODEL COMPARISON ===\n")
cat("OLS R² =", round(summary(model_ols)$r.squared, 3), "\n")
cat("OLS 95% CI: [", round(confint(model_ols)[2,1], 3), ", ", 
    round(confint(model_ols)[2,2], 3), "]\n")
cat("Robust 95% CI: [", round(confint(model_robust)[2,1], 3), ", ", 
    round(confint(model_robust)[2,2], 3), "]\n")

# 7. DIAGNOSTIC PLOTS
# Residuals vs Fitted
p1 <- ggplot(fortify(model_ols), aes(x = .fitted, y = .resid)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
  theme_minimal()

# Q-Q plot
p2 <- ggplot(fortify(model_ols), aes(sample = .stdresid)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

# Regression plot
p3 <- ggplot(rou_rep, aes(x = inspection_time, y = anti)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.3) +
  geom_smooth(method = "lmrob", se = FALSE, color = "red", linetype = "dashed") +
  labs(x = "Inspection Time", y = "Antisaccade Task Score", 
       title = "OLS (blue) vs Robust (red dashed)") +
  theme_minimal()

# Combine plots
combined_plots <- (p1 + p2) / p3 + plot_annotation(tag_levels = 'A')
print(combined_plots)

# Save plots
ggsave("diagnostic_plots.png", combined_plots, width = 12, height = 8, dpi = 300)
ggsave("regression_comparison.png", p3, width = 8, height = 6, dpi = 300)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Plots saved: diagnostic_plots.png, regression_comparison.png\n")
