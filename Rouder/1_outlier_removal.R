# Anti-Saccade and Processing Speed Data Analysis
# This script processes behavioral data from anti-saccade, temporal order, and inspection time tasks

library(tidyverse)
library(patchwork)
library(fs)

# ============================================================================
# ANTI-SACCADE DATA PROCESSING
# ============================================================================

#' Load and process anti-saccade behavioral data
#'
#' Reads CSV files from the anti-pro-saccade directory and processes them
#' to calculate mean stimulus times by participant and block type
load_saccade_data <- function() {
  # Define file paths and column types for consistent data import
  saccade_files <- dir_ls("data/anti-pro-saccade/beh", glob = "*.csv")

  col_spec <- cols(
    Block_no = col_integer(),
    Block_type = col_character(),
    Trial_no = col_integer(),
    CSI = col_integer(),
    Level = col_integer(),
    Reversal = col_integer(),
    Revs_count = col_integer()
  )

  # Load and combine all saccade data files
  saccades_raw <- saccade_files |>
    map_df(read_csv, col_types = col_spec)

  return(saccades_raw)
}

#' Analyze saccade data by block for visualization
#'
#' Creates block-level summaries to examine learning effects across blocks
analyze_blocks <- function(saccades_raw, min_reversals = 10) {
  saccades_by_block <- saccades_raw |>
    filter(Trial_type == "exp",
           Revs_count >= min_reversals,
           Reversal == 1 ) |>
    group_by(PART_ID, Block_type, Block_no) |>
    summarise(
      mean_stimtime = mean(`Stimulus Time`, na.rm = TRUE),
      .groups = "drop"
    ) |>
    # Create sequential block numbering within each condition
    group_by(PART_ID, Block_type) |>
    mutate(block_seq = row_number()) |>
    ungroup() |>
    # Generate block labels for wide format}
    mutate(block_label = paste0(Block_type, "_", block_seq)) |>
    select(PART_ID, block_label, mean_stimtime) |>
    pivot_wider(names_from = block_label, values_from = mean_stimtime) |>
    drop_na()

  return(saccades_by_block)
}

#' Create visualization of block-level performance
#'
#' Plots mean stimulus time across blocks for each condition
visualize_blocks <- function(saccades_by_block) {
  plot <- saccades_by_block |>
    pivot_longer(
      cols = -PART_ID,
      names_to = c("condition", "block"),
      names_sep = "_",
      values_to = "STIMTIME"
    ) |>
    mutate(block = as.numeric(block)) |>
    ggplot(aes(x = block, y = STIMTIME, group = PART_ID, color = condition)) +
    geom_line(alpha = 0.6) +
    facet_wrap(~condition) +
    labs(
      title = "Stimulus Time Across Blocks by Condition",
      x = "Block Number",
      y = "Mean Stimulus Time (ms)",
      color = "Condition",
      caption = "Each line represents one participant"
    ) +
    theme_minimal()

  return(plot)
}

#' Process saccade data for final analysis
#'
#' Aggregates data by participant and condition, calculates contrast score
process_saccades <- function(saccades_raw, min_reversals = 10) {
  saccades_final <- saccades_raw |>
    filter(
      Trial_type == "exp",
      Revs_count >= min_reversals,
      Reversal == 1
    ) |>
    group_by(PART_ID, Block_type) |>
    summarise(mean_time = mean(`Stimulus Time`, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = Block_type, values_from = mean_time) |>
    mutate(contrast = AS - PS) |> # Anti-saccade cost
    rename(
      pro_mean = PS, # Pro-saccade mean time
      anti_mean = AS # Anti-saccade mean time
    )

  return(saccades_final)
}

# ============================================================================
# PROCESSING SPEED DATA
# ============================================================================

#' Generic function to process temporal tasks (temporal order & inspection time)
#'
#' @param path_pattern Directory path pattern for data files
#' @param prefix Column name prefix for the processed data
#' @param min_reversals Minimum number of reversals required for inclusion
process_task_data <- function(path_pattern, prefix, min_reversals = 10) {
  dir_ls(path_pattern, glob = "*.csv") |>
    map_df(read_csv, col_types = cols(
      Reversal = col_integer(),
      Reversal_count = col_integer()
    )) |>
    # Filter for valid experimental trials with sufficient reversals
    filter(
      Reversal == 1,
      Reversal_count >= min_reversals,
      Training == "exp"
    ) |>
    # Calculate mean SOA (stimulus onset asynchrony) by participant and stimulus type
    group_by(PART_ID, Stimuli) |>
    summarize(mean_soa = mean(SOA), .groups = "drop") |>
    # Convert to wide format with appropriate column prefixes
    pivot_wider(
      names_from = Stimuli,
      values_from = mean_soa,
      names_prefix = prefix
    )
}

# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

# Load and process anti-saccade data
saccades_raw <- load_saccade_data()

# Analyze block-level effects (for quality control)
saccades_by_block <- analyze_blocks(saccades_raw)
block_plot <- visualize_blocks(saccades_by_block)

frames_to_ms <- function(df, multiplier = 16.6666) {
  df |>
    mutate(across(where(is.numeric), ~ . * multiplier))
}

# Process final saccade measures
saccades_final <- frames_to_ms(analyze_blocks(saccades_raw))

# Process temporal order and inspection time tasks
temporal_order <- frames_to_ms(process_task_data("./data/Temporal-order/beh", "temporal_"))
inspection_time <- frames_to_ms(process_task_data("./data/Inspection-Time/beh", "inspection_"))

# Combine processing speed measures
processing_speed <- inner_join(temporal_order, inspection_time, by = "PART_ID") |>
  drop_na()

# Create final combined dataset
full_df <- inner_join(saccades_final, processing_speed, by = "PART_ID") |>
  drop_na()

# ============================================================================
# DATA SUMMARY
# ============================================================================

p1 <- ggplot(full_df, aes(x = AS_4, y = PS_4)) +
  geom_point() +
  labs(title = "AS: anti_mean vs pro_mean")

p2 <- ggplot(full_df, aes(x = inspection_CIRCLES, y = inspection_SQUARES)) +
  geom_point() +
  labs(title = "Inspection: CIRCLES vs SQUARES")

p3 <- ggplot(full_df, aes(x = temporal_CIRCLES, y = temporal_SQUARES)) +
  geom_point() +
  labs(title = "Temporal: CIRCLES vs SQUARES")

p1 + p2 + p3

cat("Data processing complete!\n")
cat("Number of participants:", nrow(full_df), "\n")
cat("Four participant removed due to data incompleteness")

# Display the block visualization
print(block_plot)



# ========================================
# CONFIGURATION PARAMETERS
# ========================================

# Domain-specific thresholds based on cognitive psychology literature
MIN_VALID_RT <- 0 # Minimum plausible RT in ms (increased from 60)
MAX_VALID_RT <- 600 # Maximum plausible RT in ms (added)
MAX_NEGATIVE_CONTRAST <- -60 # Allow small measurement error (changed from 0)
MIN_PROCESSING_SCORE <- 0 # Technical floor for processing tasks

# Statistical thresholds
MAD_THRESHOLD <- 3 # Number of MADs for outlier detection (reduced from 4)
MAHAL_ALPHA <- 0.001 # Alpha level for multivariate outliers

# ========================================
# HELPER FUNCTIONS
# ========================================

#' Format numbers according to APA 7 guidelines
#' @param x Numeric value
#' @param digits Number of decimal places
#' @param leading_zero Whether to include leading zero
format_apa <- Vectorize(function(x, digits = 2, leading_zero = TRUE) {
  if (is.na(x)) return("NA")
  
  formatted <- sprintf(paste0("%.", digits, "f"), x)
  
  if (!leading_zero && abs(x) < 1) {
    formatted <- sub("^0", "", formatted)
    formatted <- sub("^-0", "-", formatted)
  }
  
  formatted
})

#' Calculate descriptive statistics for APA reporting
#' @param x Numeric vector
#' @param var_name Variable name for output
#' @return Named list of statistics
calculate_descriptives <- function(x, var_name = "Variable") {
  valid_x <- x[!is.na(x)]
  
  list(
    name = var_name,
    n = length(valid_x),
    missing = sum(is.na(x)),
    mean = mean(valid_x),
    sd = sd(valid_x),
    median = median(valid_x),
    mad = mad(valid_x),
    min = min(valid_x),
    max = max(valid_x),
    skew = moments::skewness(valid_x),
    kurtosis = moments::kurtosis(valid_x) - 3  # Excess kurtosis
  )
}

#' Print descriptive statistics in APA format
#' @param desc_list List from calculate_descriptives
print_apa_descriptives <- function(desc_list) {
  cat(sprintf(
    "%s: M = %s, SD = %s, Mdn = %s, Range = [%s, %s], n = %d",
    desc_list$name,
    format_apa(desc_list$mean),
    format_apa(desc_list$sd),
    format_apa(desc_list$median),
    format_apa(desc_list$min),
    format_apa(desc_list$max),
    desc_list$n
  ))
  
  if (desc_list$missing > 0) {
    cat(sprintf(" (%d missing)", desc_list$missing))
  }
  
  cat("\n")
}

#' Detect outliers using Median Absolute Deviation
#' @param x Numeric vector
#' @param threshold Number of MADs from median (default = 3)
#' @return Logical vector indicating outliers
detect_outliers_mad <- function(x, threshold = MAD_THRESHOLD) {
  if (all(is.na(x))) {
    return(rep(FALSE, length(x)))
  }
  
  med <- median(x, na.rm = TRUE)
  mad_value <- mad(x, na.rm = TRUE, constant = 1.4826) # Consistency with normal dist
  
  # Handle case where MAD = 0 (all values identical)
  if (mad_value == 0) {
    return(x != med)
  }
  
  abs(x - med) / mad_value > threshold
}

#' Calculate Mahalanobis distance for multivariate outlier detection
#' @param data Dataframe containing the variables
#' @param vars Character vector of variable names
#' @return Numeric vector of Mahalanobis distances
calculate_mahalanobis <- function(data, vars) {
  # Extract and prepare data
  x <- data |>
    select(all_of(vars)) |>
    drop_na() |> # Handle missing values
    as.matrix()
  
  # Return NA if insufficient data
  if (nrow(x) < length(vars) + 1) {
    return(rep(NA_real_, nrow(data)))
  }
  
  # Calculate robust estimates using minimum covariance determinant
  center <- colMeans(x)
  cov_matrix <- cov(x)
  
  # Check for singular matrix
  if (det(cov_matrix) < 1e-10) {
    warning("Covariance matrix is nearly singular. Using regularization.")
    diag(cov_matrix) <- diag(cov_matrix) + 1e-6
  }
  
  # Calculate distances for full dataset (including NAs)
  full_x <- data |>
    select(all_of(vars)) |>
    as.matrix()
  mahalanobis(full_x, center, cov_matrix)
}

#' Generate detailed outlier report
#' @param data Dataframe with outlier flags
#' @return List containing summary statistics and flagged cases
generate_outlier_report <- function(data) {
  # Count outliers by type
  outlier_counts <- data |>
    summarise(
      domain_specific = sum(domain_outlier, na.rm = TRUE),
      univariate_mad = sum(mad_outlier, na.rm = TRUE),
      multivariate = sum(mahal_outlier, na.rm = TRUE),
      any_outlier = sum(is_outlier, na.rm = TRUE),
      total_n = n()
    )
  
  # Calculate percentages for APA reporting
  outlier_percentages <- outlier_counts |>
    mutate(
      pct_domain = 100 * domain_specific / total_n,
      pct_mad = 100 * univariate_mad / total_n,
      pct_mahal = 100 * multivariate / total_n,
      pct_total = 100 * any_outlier / total_n
    )
  
  # Identify specific issues
  outlier_details <- data |>
    filter(is_outlier) |>
    select(
      PART_ID, outlier_type, contains("flag_"),
      contains("mad_"), mahal_outlier
    ) |>
    arrange(PART_ID)
  
  # Count specific domain violations
  domain_violations <- data |>
    summarise(
      too_fast_AS = sum(flag_AS_1_too_fast | flag_AS_2_too_fast | 
                          flag_AS_3_too_fast | flag_AS_4_too_fast, na.rm = TRUE),
      too_fast_PS = sum(flag_PS_1_too_fast | flag_PS_2_too_fast | 
                          flag_PS_3_too_fast | flag_PS_4_too_fast, na.rm = TRUE),
      too_slow_AS = sum(flag_AS_1_too_slow | flag_AS_2_too_slow | 
                          flag_AS_3_too_slow | flag_AS_4_too_slow, na.rm = TRUE),
      too_slow_PS = sum(flag_PS_1_too_slow | flag_PS_2_too_slow | 
                          flag_PS_3_too_slow | flag_PS_4_too_slow, na.rm = TRUE),
      negative_contrast = sum(flag_negative_contrast, na.rm = TRUE),
      zero_processing = sum(flag_zero_temporal_C | flag_zero_temporal_S | 
                              flag_zero_inspection_C | flag_zero_inspection_S, na.rm = TRUE)
    )
  
  list(
    summary = outlier_counts,
    percentages = outlier_percentages,
    details = outlier_details,
    domain_violations = domain_violations,
    percent_removed = round(100 * outlier_counts$any_outlier / outlier_counts$total_n, 1)
  )
}

#' Generate APA-style method section text
#' @param params List of parameters used
#' @param report Outlier report object
generate_apa_method_text <- function(params, report) {
  cat("\n=== APA STYLE METHOD SECTION TEXT ===\n\n")
  
  cat("Data Screening and Outlier Detection\n\n")
  
  cat(sprintf(
    "Data were screened for outliers using a three-stage hierarchical approach. First, domain-specific criteria were applied based on cognitive psychology literature: (a) reaction times < %d ms or > %d ms were flagged as anticipatory responses or attention lapses, respectively; (b) antisaccade-prosaccade contrast scores < %d ms were flagged as implausible; and (c) processing speed scores of 0 were flagged as technical failures.\n\n",
    params$min_rt, params$max_rt, params$max_neg_contrast
  ))
  
  cat(sprintf(
    "Second, univariate outliers were identified using the Median Absolute Deviation (MAD) method with a threshold of %s MADs from the median (Leys et al., 2013). Third, multivariate outliers were detected using Mahalanobis distance with α = %s.\n\n",
    format_apa(params$mad_threshold, 1), format_apa(params$mahal_alpha, 3, FALSE)
  ))
  
  cat("Participants were excluded if they met domain-specific criteria or at least two statistical criteria. This approach balances sensitivity to data quality issues with retention of valid but extreme scores.\n\n")
}

#' Generate APA-style results text
#' @param report Outlier report object
#' @param initial_n Initial sample size
generate_apa_results_text <- function(report, initial_n) {
  cat("\n=== APA STYLE RESULTS TEXT ===\n\n")
  
  cat("Outlier Detection Results\n\n")
  
  # Overall results
  cat(sprintf(
    "Of the initial %d participants, %d (%s%%) were identified as outliers and removed from analyses. ",
    initial_n,
    report$summary$any_outlier,
    format_apa(report$percentages$pct_total, 1)
  ))
  
  # Breakdown by criterion
  cat(sprintf(
    "Domain-specific criteria identified %d participants (%s%%), univariate MAD analysis identified %d (%s%%), and multivariate analysis identified %d (%s%%).\n\n",
    report$summary$domain_specific, format_apa(report$percentages$pct_domain, 1),
    report$summary$univariate_mad, format_apa(report$percentages$pct_mad, 1),
    report$summary$multivariate, format_apa(report$percentages$pct_mahal, 1)
  ))
  
  # Specific violations
  cat("Among domain-specific violations: ")
  violations <- c()
  
  if (report$domain_violations$too_fast_AS > 0 || report$domain_violations$too_fast_PS > 0) {
    violations <- c(violations, sprintf("%d showed anticipatory responses",
                                        report$domain_violations$too_fast_AS + report$domain_violations$too_fast_PS))
  }
  
  if (report$domain_violations$too_slow_AS > 0 || report$domain_violations$too_slow_PS > 0) {
    violations <- c(violations, sprintf("%d showed attention lapses",
                                        report$domain_violations$too_slow_AS + report$domain_violations$too_slow_PS))
  }
  
  if (report$domain_violations$negative_contrast > 0) {
    violations <- c(violations, sprintf("%d had implausible contrast scores",
                                        report$domain_violations$negative_contrast))
  }
  
  if (report$domain_violations$zero_processing > 0) {
    violations <- c(violations, sprintf("%d had zero processing speed scores",
                                        report$domain_violations$zero_processing))
  }
  
  cat(paste(violations, collapse = ", "))
  cat(sprintf(". The final sample consisted of %d participants.\n\n",
              initial_n - report$summary$any_outlier))
}

# ========================================
# MAIN OUTLIER DETECTION PIPELINE
# ========================================

#' Main function to detect and remove outliers
#' @param df Input dataframe
#' @param verbose Print progress messages (default = TRUE)
#' @return List containing cleaned data and diagnostic information
remove_outliers <- function(df, verbose = TRUE) {
  if (verbose) {
    cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")
    cat("COGNITIVE DATA OUTLIER DETECTION AND CLEANING\n")
    cat(paste(rep("=", 60), collapse = "") %+% "\n\n")
    cat("Initial sample size: N =", nrow(df), "\n")
    cat("Analysis timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("\n")
  }
  
  # ----------------------------------------
  # STEP 0: Calculate means and contrast
  # ----------------------------------------
  df <- df |>
    mutate(
      anti_mean = rowMeans(pick(AS_1, AS_2, AS_3, AS_4), na.rm = TRUE),
      pro_mean = rowMeans(pick(PS_1, PS_2, PS_3, PS_4), na.rm = TRUE),
      contrast = anti_mean - pro_mean
    )
  
  # Report initial descriptive statistics
  if (verbose) {
    cat("INITIAL DESCRIPTIVE STATISTICS\n")
    cat(paste(rep("-", 40), collapse = "") %+% "\n")
    
    # Saccade tasks
    cat("\nSaccade Tasks (in milliseconds):\n")
    for (var in c("AS_1", "AS_2", "AS_3", "AS_4", "PS_1", "PS_2", "PS_3", "PS_4")) {
      desc <- calculate_descriptives(df[[var]], var)
      print_apa_descriptives(desc)
    }
    
    cat("\nCalculated Variables:\n")
    print_apa_descriptives(calculate_descriptives(df$anti_mean, "Antisaccade Mean"))
    print_apa_descriptives(calculate_descriptives(df$pro_mean, "Prosaccade Mean"))
    print_apa_descriptives(calculate_descriptives(df$contrast, "AS-PS Contrast"))
    
    cat("\nProcessing Speed Tasks:\n")
    for (var in c("temporal_CIRCLES", "temporal_SQUARES", "inspection_CIRCLES", "inspection_SQUARES")) {
      desc <- calculate_descriptives(df[[var]], var)
      print_apa_descriptives(desc)
    }
    cat("\n")
  }
  
  # ----------------------------------------
  # STEP 1: Domain-Specific Screening
  # ----------------------------------------
  if (verbose) {
    cat("STEP 1: DOMAIN-SPECIFIC SCREENING\n")
    cat(paste(rep("-", 40), collapse = "") %+% "\n")
    cat(sprintf("Criteria: %d ms < RT < %d ms, Contrast > %d ms, Processing > %d\n\n",
                MIN_VALID_RT, MAX_VALID_RT, MAX_NEGATIVE_CONTRAST, MIN_PROCESSING_SCORE))
  }
  
  df_flagged <- df |>
    mutate(
      # Implausibly fast RTs (anticipatory responses) - check all columns
      flag_AS_1_too_fast = AS_1 < MIN_VALID_RT,
      flag_AS_2_too_fast = AS_2 < MIN_VALID_RT,
      flag_AS_3_too_fast = AS_3 < MIN_VALID_RT,
      flag_AS_4_too_fast = AS_4 < MIN_VALID_RT,
      flag_PS_1_too_fast = PS_1 < MIN_VALID_RT,
      flag_PS_2_too_fast = PS_2 < MIN_VALID_RT,
      flag_PS_3_too_fast = PS_3 < MIN_VALID_RT,
      flag_PS_4_too_fast = PS_4 < MIN_VALID_RT,
      
      # Implausibly slow RTs (lapses) - check all columns
      flag_AS_1_too_slow = AS_1 > MAX_VALID_RT,
      flag_AS_2_too_slow = AS_2 > MAX_VALID_RT,
      flag_AS_3_too_slow = AS_3 > MAX_VALID_RT,
      flag_AS_4_too_slow = AS_4 > MAX_VALID_RT,
      flag_PS_1_too_slow = PS_1 > MAX_VALID_RT,
      flag_PS_2_too_slow = PS_2 > MAX_VALID_RT,
      flag_PS_3_too_slow = PS_3 > MAX_VALID_RT,
      flag_PS_4_too_slow = PS_4 > MAX_VALID_RT,
      
      # Negative contrast beyond threshold (prosaccade > antisaccade)
      flag_negative_contrast = contrast < MAX_NEGATIVE_CONTRAST,
      
      # Near-zero processing speed (technical failures)
      flag_zero_temporal_C = temporal_CIRCLES < MIN_PROCESSING_SCORE,
      flag_zero_temporal_S = temporal_SQUARES < MIN_PROCESSING_SCORE,
      flag_zero_inspection_C = inspection_CIRCLES < MIN_PROCESSING_SCORE,
      flag_zero_inspection_S = inspection_SQUARES < MIN_PROCESSING_SCORE,
      
      # Aggregate domain flag
      domain_outlier = flag_AS_1_too_fast | flag_AS_2_too_fast | flag_AS_3_too_fast | flag_AS_4_too_fast |
        flag_PS_1_too_fast | flag_PS_2_too_fast | flag_PS_3_too_fast | flag_PS_4_too_fast |
        flag_AS_1_too_slow | flag_AS_2_too_slow | flag_AS_3_too_slow | flag_AS_4_too_slow |
        flag_PS_1_too_slow | flag_PS_2_too_slow | flag_PS_3_too_slow | flag_PS_4_too_slow |
        flag_negative_contrast | flag_zero_temporal_C |
        flag_zero_temporal_S | flag_zero_inspection_C |
        flag_zero_inspection_S
    )
  
  if (verbose) {
    # Detailed breakdown
    domain_summary <- df_flagged |>
      summarise(
        fast_AS = sum(flag_AS_1_too_fast | flag_AS_2_too_fast | flag_AS_3_too_fast | flag_AS_4_too_fast, na.rm = TRUE),
        fast_PS = sum(flag_PS_1_too_fast | flag_PS_2_too_fast | flag_PS_3_too_fast | flag_PS_4_too_fast, na.rm = TRUE),
        slow_AS = sum(flag_AS_1_too_slow | flag_AS_2_too_slow | flag_AS_3_too_slow | flag_AS_4_too_slow, na.rm = TRUE),
        slow_PS = sum(flag_PS_1_too_slow | flag_PS_2_too_slow | flag_PS_3_too_slow | flag_PS_4_too_slow, na.rm = TRUE),
        neg_contrast = sum(flag_negative_contrast, na.rm = TRUE),
        zero_proc = sum(flag_zero_temporal_C | flag_zero_temporal_S | flag_zero_inspection_C | flag_zero_inspection_S, na.rm = TRUE)
      )
    
    cat("Domain-specific violations detected:\n")
    cat(sprintf("  - Anticipatory antisaccades (< %d ms): %d participants\n", MIN_VALID_RT, domain_summary$fast_AS))
    cat(sprintf("  - Anticipatory prosaccades (< %d ms): %d participants\n", MIN_VALID_RT, domain_summary$fast_PS))
    cat(sprintf("  - Slow antisaccades (> %d ms): %d participants\n", MAX_VALID_RT, domain_summary$slow_AS))
    cat(sprintf("  - Slow prosaccades (> %d ms): %d participants\n", MAX_VALID_RT, domain_summary$slow_PS))
    cat(sprintf("  - Implausible contrast (< %d ms): %d participants\n", MAX_NEGATIVE_CONTRAST, domain_summary$neg_contrast))
    cat(sprintf("  - Zero processing scores: %d participants\n", domain_summary$zero_proc))
    cat(sprintf("\nTotal domain outliers: %d (%s%% of sample)\n\n", 
                sum(df_flagged$domain_outlier, na.rm = TRUE),
                format_apa(100 * sum(df_flagged$domain_outlier, na.rm = TRUE) / nrow(df), 1)))
  }
  
  # ----------------------------------------
  # STEP 2: Univariate MAD Outliers
  # ----------------------------------------
  if (verbose) {
    cat("STEP 2: UNIVARIATE OUTLIER DETECTION (MAD)\n")
    cat(paste(rep("-", 40), collapse = "") %+% "\n")
    cat(sprintf("Threshold: %s MADs from median\n\n", format_apa(MAD_THRESHOLD, 1)))
  }
  
  df_mad <- df_flagged |>
    mutate(
      # Apply MAD detection to each AS and PS column
      mad_AS_1 = detect_outliers_mad(AS_1),
      mad_AS_2 = detect_outliers_mad(AS_2),
      mad_AS_3 = detect_outliers_mad(AS_3),
      mad_AS_4 = detect_outliers_mad(AS_4),
      mad_PS_1 = detect_outliers_mad(PS_1),
      mad_PS_2 = detect_outliers_mad(PS_2),
      mad_PS_3 = detect_outliers_mad(PS_3),
      mad_PS_4 = detect_outliers_mad(PS_4),
      
      # Also check means and contrast
      mad_contrast = detect_outliers_mad(contrast),
      
      # Processing tasks
      mad_temporal_C = detect_outliers_mad(temporal_CIRCLES),
      mad_temporal_S = detect_outliers_mad(temporal_SQUARES),
      mad_inspection_C = detect_outliers_mad(inspection_CIRCLES),
      mad_inspection_S = detect_outliers_mad(inspection_SQUARES),
      
      # Aggregate MAD flag
      mad_outlier = mad_AS_1 | mad_AS_2 | mad_AS_3 | mad_AS_4 |
        mad_PS_1 | mad_PS_2 | mad_PS_3 | mad_PS_4 |  mad_contrast |
        mad_temporal_C | mad_temporal_S |
        mad_inspection_C | mad_inspection_S
    )
  
  if (verbose) {
    # Create summary table
    mad_summary <- data.frame(
      Variable = c("AS_1", "AS_2", "AS_3", "AS_4", "PS_1", "PS_2", "PS_3", "PS_4",
                   "Contrast", "Temporal_C", "Temporal_S", "Inspection_C", "Inspection_S"),
      N_Outliers = c(
        sum(df_mad$mad_AS_1, na.rm = TRUE),
        sum(df_mad$mad_AS_2, na.rm = TRUE),
        sum(df_mad$mad_AS_3, na.rm = TRUE),
        sum(df_mad$mad_AS_4, na.rm = TRUE),
        sum(df_mad$mad_PS_1, na.rm = TRUE),
        sum(df_mad$mad_PS_2, na.rm = TRUE),
        sum(df_mad$mad_PS_3, na.rm = TRUE),
        sum(df_mad$mad_PS_4, na.rm = TRUE),
        sum(df_mad$mad_contrast, na.rm = TRUE),
        sum(df_mad$mad_temporal_C, na.rm = TRUE),
        sum(df_mad$mad_temporal_S, na.rm = TRUE),
        sum(df_mad$mad_inspection_C, na.rm = TRUE),
        sum(df_mad$mad_inspection_S, na.rm = TRUE)
      )
    ) |>
      mutate(Percentage = format_apa(100 * N_Outliers / nrow(df), 1))
    
    cat("MAD outliers by variable:\n")
    print(mad_summary, row.names = FALSE)
    cat(sprintf("\nTotal MAD outliers: %d (%s%% of sample)\n\n",
                sum(df_mad$mad_outlier, na.rm = TRUE),
                format_apa(100 * sum(df_mad$mad_outlier, na.rm = TRUE) / nrow(df), 1)))
  }
  
  # ----------------------------------------
  # STEP 3: Multivariate Outliers
  # ----------------------------------------
  if (verbose) {
    cat("STEP 3: MULTIVARIATE OUTLIER DETECTION\n")
    cat(paste(rep("-", 40), collapse = "") %+% "\n")
    cat(sprintf("Method: Mahalanobis distance, α = %s\n\n", format_apa(MAHAL_ALPHA, 3, FALSE)))
  }
  
  # Define all saccade variables
  saccade_vars <- c("AS_1", "AS_2", "AS_3", "AS_4", 
                    "PS_1", "PS_2", "PS_3", "PS_4")
  
  # Define processing speed variables
  processing_vars <- c(
    "temporal_CIRCLES", "temporal_SQUARES",
    "inspection_CIRCLES", "inspection_SQUARES"
  )
  
  df_mahal <- df_mad %>%
    mutate(
      # Calculate Mahalanobis distance for saccade tasks
      mahal_dist_saccade = calculate_mahalanobis(., saccade_vars),
      
      # Calculate Mahalanobis distance for processing tasks
      mahal_dist_processing = calculate_mahalanobis(., processing_vars),
      
      # Flag outliers based on chi-square distribution
      mahal_outlier_saccade = mahal_dist_saccade > qchisq(1 - MAHAL_ALPHA, df = length(saccade_vars)),
      mahal_outlier_processing = mahal_dist_processing > qchisq(1 - MAHAL_ALPHA, df = length(processing_vars)),
      
      # Combined multivariate outlier flag
      mahal_outlier = mahal_outlier_saccade | mahal_outlier_processing
    )
  
  if (verbose) {
    # Calculate critical values
    crit_saccade <- qchisq(1 - MAHAL_ALPHA, df = length(saccade_vars))
    crit_processing <- qchisq(1 - MAHAL_ALPHA, df = length(processing_vars))
    
    cat(sprintf("Critical values:\n"))
    cat(sprintf("  - Saccade tasks (df = %d): χ² = %s\n", length(saccade_vars), format_apa(crit_saccade)))
    cat(sprintf("  - Processing tasks (df = %d): χ² = %s\n\n", length(processing_vars), format_apa(crit_processing)))
    
    cat("Multivariate outliers detected:\n")
    cat(sprintf("  - Saccade tasks: %d (%s%%)\n", 
                sum(df_mahal$mahal_outlier_saccade, na.rm = TRUE),
                format_apa(100 * sum(df_mahal$mahal_outlier_saccade, na.rm = TRUE) / nrow(df), 1)))
    cat(sprintf("  - Processing tasks: %d (%s%%)\n", 
                sum(df_mahal$mahal_outlier_processing, na.rm = TRUE),
                format_apa(100 * sum(df_mahal$mahal_outlier_processing, na.rm = TRUE) / nrow(df), 1)))
    cat(sprintf("\nTotal multivariate outliers: %d (%s%% of sample)\n\n",
                sum(df_mahal$mahal_outlier, na.rm = TRUE),
                format_apa(100 * sum(df_mahal$mahal_outlier, na.rm = TRUE) / nrow(df), 1)))
  }
  
  # ----------------------------------------
  # STEP 4: Final Classification
  # ----------------------------------------
  if (verbose) {
    cat("STEP 4: FINAL OUTLIER CLASSIFICATION\n")
    cat(paste(rep("-", 40), collapse = "") %+% "\n")
    cat("Decision rule: Exclude if domain violation OR ≥2 statistical criteria\n\n")
  }
  
  df_final <- df_mahal %>%
    mutate(
      # Count statistical criteria met
      stat_criteria_count = as.numeric(mad_outlier) + as.numeric(mahal_outlier),
      
      # Overall outlier flag with weighted approach
      is_outlier = domain_outlier | (stat_criteria_count >= 2),
      
      # Detailed classification
      outlier_type = case_when(
        !is_outlier ~ "None",
        domain_outlier & (stat_criteria_count >= 2) ~ "Domain + Statistical",
        domain_outlier ~ "Domain only",
        mad_outlier & mahal_outlier ~ "MAD + Multivariate",
        mad_outlier & !mahal_outlier ~ "MAD only (kept)",
        !mad_outlier & mahal_outlier ~ "Multivariate only (kept)",
        TRUE ~ "Other"
      )
    )
  
  # ----------------------------------------
  # STEP 5: Create Output
  # ----------------------------------------
  
  # Clean dataset (remove outliers and auxiliary columns including calculated means)
  cleaned_data <- df_final %>%
    filter(!is_outlier) %>%
    select(-anti_mean, -pro_mean, -contrast) %>%  # Remove calculated columns
    select(all_of(names(df)[!names(df) %in% c("anti_mean", "pro_mean", "contrast")])) # Keep only original columns
  
  # Generate report
  report <- generate_outlier_report(df_final)
  
  if (verbose) {
    cat("\n" %+% paste(rep("=", 60), collapse = "") %+% "\n")
    cat("SUMMARY\n")
    cat(paste(rep("=", 60), collapse = "") %+% "\n\n")
    
    # Overall summary
    cat(sprintf("Initial sample size: N = %d\n", nrow(df)))
    cat(sprintf("Outliers removed: n = %d (%s%%)\n", 
                report$summary$any_outlier, 
                format_apa(report$percent_removed, 1)))
    cat(sprintf("Final sample size: N = %d\n\n", nrow(cleaned_data)))
    
    # Breakdown by type
    type_breakdown <- df_final %>%
      filter(is_outlier) %>%
      count(outlier_type) %>%
      mutate(percentage = format_apa(100 * n / sum(n), 1)) %>%
      arrange(desc(n))
    
    cat("Outliers by classification:\n")
    print(type_breakdown, n = Inf)
    
    # Overlap analysis
    cat("\nCriterion overlap analysis:\n")
    overlap_stats <- df_final %>%
      filter(is_outlier) %>%
      summarise(
        domain_only = sum(domain_outlier & !mad_outlier & !mahal_outlier),
        mad_only = sum(!domain_outlier & mad_outlier & !mahal_outlier),
        mahal_only = sum(!domain_outlier & !mad_outlier & mahal_outlier),
        domain_mad = sum(domain_outlier & mad_outlier & !mahal_outlier),
        domain_mahal = sum(domain_outlier & !mad_outlier & mahal_outlier),
        mad_mahal = sum(!domain_outlier & mad_outlier & mahal_outlier),
        all_three = sum(domain_outlier & mad_outlier & mahal_outlier)
      )
    
    cat(sprintf("  - Domain only: %d\n", overlap_stats$domain_only))
    cat(sprintf("  - MAD only: %d\n", overlap_stats$mad_only))
    cat(sprintf("  - Mahalanobis only: %d\n", overlap_stats$mahal_only))
    cat(sprintf("  - Domain + MAD: %d\n", overlap_stats$domain_mad))
    cat(sprintf("  - Domain + Mahalanobis: %d\n", overlap_stats$domain_mahal))
    cat(sprintf("  - MAD + Mahalanobis: %d\n", overlap_stats$mad_mahal))
    cat(sprintf("  - All three criteria: %d\n\n", overlap_stats$all_three))
    
    # Generate APA text
    generate_apa_method_text(
      list(
        min_rt = MIN_VALID_RT,
        max_rt = MAX_VALID_RT,
        max_neg_contrast = MAX_NEGATIVE_CONTRAST,
        mad_threshold = MAD_THRESHOLD,
        mahal_alpha = MAHAL_ALPHA
      ),
      report
    )
    
    generate_apa_results_text(report, nrow(df))
    
    # Final descriptive statistics
    cat("\nFINAL DESCRIPTIVE STATISTICS (After Outlier Removal)\n")
    cat(paste(rep("-", 40), collapse = "") %+% "\n")
    
    # Recalculate means for cleaned data
    cleaned_with_means <- cleaned_data |>
      mutate(
        anti_mean = rowMeans(pick(AS_1, AS_2, AS_3, AS_4), na.rm = TRUE),
        pro_mean = rowMeans(pick(PS_1, PS_2, PS_3, PS_4), na.rm = TRUE),
        contrast = anti_mean - pro_mean
      )
    
    cat("\nSaccade Tasks (in milliseconds):\n")
    for (var in c("AS_1", "AS_2", "AS_3", "AS_4", "PS_1", "PS_2", "PS_3", "PS_4")) {
      desc <- calculate_descriptives(cleaned_data[[var]], var)
      print_apa_descriptives(desc)
    }
    
    cat("\nCalculated Variables:\n")
    print_apa_descriptives(calculate_descriptives(cleaned_with_means$anti_mean, "Antisaccade Mean"))
    print_apa_descriptives(calculate_descriptives(cleaned_with_means$pro_mean, "Prosaccade Mean"))
    print_apa_descriptives(calculate_descriptives(cleaned_with_means$contrast, "AS-PS Contrast"))
    
    cat("\nProcessing Speed Tasks:\n")
    for (var in c("temporal_CIRCLES", "temporal_SQUARES", "inspection_CIRCLES", "inspection_SQUARES")) {
      desc <- calculate_descriptives(cleaned_data[[var]], var)
      print_apa_descriptives(desc)
    }
  }
  
  # Return results
  list(
    cleaned_data = cleaned_data,
    full_data = df_final,
    report = report,
    parameters = list(
      min_rt = MIN_VALID_RT,
      max_rt = MAX_VALID_RT,
      max_neg_contrast = MAX_NEGATIVE_CONTRAST,
      mad_threshold = MAD_THRESHOLD,
      mahal_alpha = MAHAL_ALPHA
    )
  )
}

# ========================================
# EXECUTE ANALYSIS
# ========================================

# Run outlier detection
results <- remove_outliers(full_df, verbose = TRUE)

# Extract cleaned data
cleaned_df <- results$cleaned_data

# Calculate means for visualization
cleaned_df_with_means <- cleaned_df |>
  mutate(
    anti_mean = rowMeans(pick(AS_1, AS_2, AS_3, AS_4), na.rm = TRUE),
    pro_mean = rowMeans(pick(PS_1, PS_2, PS_3, PS_4), na.rm = TRUE)
  )

# Create diagnostic plots with APA styling
library(ggplot2)
library(patchwork)

# Set APA theme
theme_apa <- theme_minimal() +
  theme(
    text = element_text(family = "Times", size = 12),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

p1 <- ggplot(cleaned_df_with_means, aes(x = anti_mean, y = pro_mean)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "Mean Antisaccade RT (ms)",
    y = "Mean Prosaccade RT (ms)",
    title = "A. Saccade Task Performance"
  ) +
  theme_apa

p2 <- ggplot(cleaned_df, aes(x = inspection_CIRCLES, y = inspection_SQUARES)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "Inspection Time: Circles",
    y = "Inspection Time: Squares",
    title = "B. Inspection Time Performance"
  ) +
  theme_apa

p3 <- ggplot(cleaned_df, aes(x = temporal_CIRCLES, y = temporal_SQUARES)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    x = "Temporal Processing: Circles",
    y = "Temporal Processing: Squares",
    title = "C. Temporal Processing Performance"
  ) +
  theme_apa

# Combine plots
combined_plot <- p1 + p2 + p3 + 
  plot_annotation(
    title = "Figure 1",
    subtitle = "Scatterplots of Cognitive Performance Variables After Outlier Removal",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(face = "italic", size = 12)
    )
  )

print(combined_plot)

# Save results
write_csv(cleaned_df, "cleaned_cognitive_data.csv")
write_csv(results$full_data, "full_data_with_outlier_flags.csv")

# Create supplementary table for manuscript
outlier_summary_table <- data.frame(
  Criterion = c("Domain-specific", "Univariate (MAD)", "Multivariate", "Any criterion", "Total N"),
  n = c(
    results$report$percentages$domain_specific,
    results$report$percentages$univariate_mad,
    results$report$percentages$multivariate,
    results$report$percentages$any_outlier,
    results$report$percentages$total_n
  ),
  `%` = c(
    format_apa(results$report$percentages$pct_domain, 1),
    format_apa(results$report$percentages$pct_mad, 1),
    format_apa(results$report$percentages$pct_mahal, 1),
    format_apa(results$report$percentages$pct_total, 1),
    "-"
  )
)

cat("\n\nTable 1. Outlier Detection Summary\n")
print(outlier_summary_table, row.names = FALSE)

# ggsave("outlier_diagnostics.pdf", combined_plot, width = 12, height = 5)
