library(testthat)
library(dplyr)
source("functions2.r")

# Create test data
create_test_data <- function(n_ids = 2, points_per_id = 5) {
  data.frame(
    id = rep(1:n_ids, each = points_per_id),
    time = rep(1:points_per_id, n_ids),
    corr = rnorm(n_ids * points_per_id)  # Random correlations
  )
}

# Main test suite
test_that("preprocess function handles basic operations correctly", {
  # Setup basic test data
  test_data <- data.frame(
    id = c(1, 1, 1, 2, 2, 2),
    time = c(1, 2, 3, 1, 2, 3),
    corr = c(0.5, 1.0, 1.5, -0.5, 0.0, 0.5)
  )
  
  # Test 1: Basic functionality without padding or outlier removal
  result <- preprocess(test_data, outlier = Inf, padlen = 0)
  expect_equal(nrow(result), nrow(test_data))
  expect_true(all(c("id", "time", "corr", "padding") %in% names(result)))
  expect_true(all(result$padding == 0))
  
  # Test 2: Normalization
  # Each group should have mean ≈ 0 and sd ≈ 1
  group_stats <- result %>%
    group_by(id) %>%
    summarise(
      mean_corr = mean(corr),
      sd_corr = sd(corr)
    )
  expect_true(all(abs(group_stats$mean_corr) < 1e-10))  # Mean should be very close to 0
  expect_true(all(abs(group_stats$sd_corr - 1) < 1e-10))  # SD should be very close to 1
})

test_that("padding is applied correctly", {
  # Setup simple test data
  test_data <- data.frame(
    id = c(1, 1, 1),
    time = c(1, 2, 3),
    corr = c(0.5, 1.0, 1.5)
  )
  
  # Test with padding
  padlen <- 2
  result <- preprocess(test_data, padlen = padlen)
  
  # Check padding length
  expect_equal(nrow(result), nrow(test_data) + 4)  # 2 padding points before and after
  
  # Check padding values
  padding_points <- result$padding == 1
  expect_equal(sum(padding_points), 4)  # Should have 4 padding points
  expect_true(all(result$corr[padding_points] == 0))  # Padding should have zero correlation
  
  # Check time sequence
  time_diffs <- diff(result$time)
  expect_true(all(abs(diff(time_diffs)) < 1e-10))  # Time steps should be uniform
})

test_that("outlier removal works correctly", {
  # Setup test data with outliers
  test_data <- data.frame(
    id = c(1, 1, 1, 1, 1),
    time = 1:5,
    corr = c(0.5, 10.0, 1.0, -8.0, 0.8)  # 10.0 and -8.0 are outliers
  )
  
  # Test outlier removal
  outlier_threshold <- 2
  result <- preprocess(test_data, outlier = outlier_threshold, padlen = 0)
  
  # Check if outliers were removed
  expect_true(all(abs(result$corr) < outlier_threshold))
  expect_equal(nrow(result), 3)  # Should have removed 2 points
})

test_that("function handles edge cases", {
  # Test 1: Empty dataframe
  expect_error(preprocess(data.frame()))
  
  # Test 2: Missing required columns
  expect_error(preprocess(data.frame(id = 1, time = 1)))
  
  # Test 3: Single observation
  single_obs <- data.frame(id = 1, time = 1, corr = 0.5)
  expect_error(preprocess(single_obs), NA)  # Should not error
  
  # Test 4: All identical values
  same_values <- data.frame(
    id = c(1, 1, 1),
    time = 1:3,
    corr = rep(1, 3)
  )
  expect_warning(preprocess(same_values))  # Should warn about zero standard deviation
  
  # Test 5: NA values
  na_data <- data.frame(
    id = c(1, 1, 1),
    time = 1:3,
    corr = c(0.5, NA, 1.0)
  )
  expect_error(preprocess(na_data))  # Should handle or error on NA values
})

test_that("sorting is correct", {
  # Setup unsorted test data
  test_data <- data.frame(
    id = c(2, 1, 2, 1, 2, 1),
    time = c(2, 3, 1, 1, 3, 2),
    corr = rnorm(6)
  )
  
  result <- preprocess(test_data, padlen = 0)
  
  # Check if sorted by id then time
  expect_true(all(diff(result$id) >= 0))  # IDs should be non-decreasing
  
  # Within each ID, time should be strictly increasing
  by_id <- split(result$time, result$id)
  expect_true(all(sapply(by_id, function(x) all(diff(x) > 0))))
})

test_that("performance with larger datasets", {
  # Skip on CRAN to avoid long tests
  skip_on_cran()
  
  # Create larger test dataset
  large_data <- create_test_data(n_ids = 100, points_per_id = 1000)
  
  # Test execution time
  start_time <- Sys.time()
  result <- preprocess(large_data, padlen = 50)
  end_time <- Sys.time()
  
  # Test should complete in reasonable time (adjust threshold as needed)
  expect_true(as.numeric(end_time - start_time) < 10)  # Should complete in under 10 seconds
  
  # Check result size
  expected_rows <- nrow(large_data) + (100 * 50 * 2)  # Original rows + padding for each ID
  expect_equal(nrow(result), expected_rows)
})

# Run all tests
test_results <- test_dir("tests", reporter = "summary")