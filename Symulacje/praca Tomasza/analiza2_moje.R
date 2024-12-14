library(ggplot2)
library(signal)
library(plyr)
library(fs)
library(dplyr)
library(reshape2)
library(tidyverse)
library(detrendeR)
library(progressr)
library(future)
library(furrr)

source("functions2_moje.R")

#repo_low <- "https://raw.githubusercontent.com/bartekkroczek/Behavioral-oscillations/master/data/low_resolution/"
#files_low <- read.table("files_low.txt")$V1
freq_low <- 1 / 60

# d1 <- NULL
#for (i in files_low)
#  d1 <- rbind(d1, read.csv(paste0(repo_low, i)))
#d1$Corr = as.numeric(as.factor(d1$Corr)) - 1

csv_files <- dir_ls(path = "../../Experiment_1/Data/", glob = "*.csv")

df_raw <- csv_files %>%
  map_df(~read_csv(., show_col_types = FALSE))


#d1 <- dplyr::rename(df, id = "PART_ID")
#d1 <- dplyr::rename(d1, corr = "Corr")
#d1 <- subset(d1, Trial_type == "experiment")
#d1atemp <- aggregate(corr ~ id, d1, FUN = mean)
#d1a <- aggregate(corr ~ CSI * id, d1, FUN = mean)
#d1a <- subset(d1a, id %in% d1atemp$id[d1atemp$corr > .6])
#d1a$time <- with(d1a, CSI * freq_low)


df <- df_raw %>%
  rename(id = PART_ID, corr = Corr) |>
  dplyr::filter(Trial_type == "experiment") |>
  group_by(id) |>
  mutate(mean_corr = mean(corr)) |>
  dplyr::filter(mean_corr > 0.6) |>
  group_by(CSI, id) |>
  summarize(corr = mean(corr), .groups = "drop") |>
  mutate(time = CSI * freq_low) |> 
  arrange(id, CSI) |>
  ungroup()

roll_mean <- (function(x, n = 3)
  stats::filter(x, rep(1 / n, n), sides = 2))

process_data <- function(df_raw, freq_low, corr_threshold = 0.6) {
  df_processed <- df_raw |>
    rename(id = PART_ID, corr = Corr) |># Rename columns for clarity
    dplyr::filter(Trial_type == "experiment") |># Keep only experimental trials
    group_by(id) |># Remove participants with low average performance
    dplyr::filter(mean(corr) > corr_threshold) |>
    group_by(CSI, id) |> # Calculate mean corr for each participant and CSI combination
    summarise(corr = mean(corr), .groups = "drop") |>
    mutate(time = CSI * freq_low) |> # Calculate time based on CSI and frequency
    arrange(id, CSI) |> # Sort data for rolling mean calculation
    group_by(id) |> # Calculate rolling mean correlation by participant
    mutate(rolling_corr = roll_mean(corr)) |>
    group_by(id) |> # Replace missing values with participant-specific means
    mutate(corr = if_else(is.na(rolling_corr), mean(rolling_corr, na.rm = TRUE),
                          rolling_corr)) |>
    select(id, CSI, time, corr) |> # Select final columns and ungroup
    ungroup()
  return(df_processed)
}

df2 <- process_data(df_raw, freq_low)



fftwhole <- function(data, FUN = mean) {
  temp <- remove_trend(
    hanning_window(preprocess(data, 4, 0), .(id)),
    type = "pol",
    degree = 2 * 2 * diff(range(data$time))
  )
  temp$len <- diff(range(data$time))
  aggregate(fft ~ freq, FUN = FUN, fourier(temp, .(id)))
}

permutevals <- function(data, columns) {
  ret <- ddply(data, columns, function(df) {
    cbind(df$time, melt(sample(df$corr)))
  })
  ret <- dplyr::rename(ret, corr = "value")
  ret <- dplyr::rename(ret, time = "df$time")
  ret
}

makedistr <- function(data, n = 10000, columns = c("id"), FUN = fftwhole) {
  freq <- data %>%
    group_by(across(all_of(columns))) %>%
    summarize(length = n()) %>%
    pull(length) %>%
    max() %>%
    seq_len() %>%
    {(. - 1) / diff(range(data$time))}
  
  # Set up progress reporting
  handlers(handler_progress(format = ":spin :current/:total (:percent) [:bar] ETA: :eta"))
  
  # Use future_map to parallelize the computation
  plan(multisession)
  
  with_progress({
    p <- progressor(steps = n)
    
    result <- future_map(1:n, function(i) {
      p()
      data %>%
        permutevals(columns) %>%
        FUN() %>%
        pull(fft)
    }, .options = furrr_options(seed = TRUE)) %>%
      bind_cols()
  })
  
  # Convert the result to a matrix and set row names
  result <- as.matrix(result)
  rownames(result) <- freq
  
  return(result)
}

plotci <- function(cimat,
                   estdf,
                   maxfreq = 31,
                   alpha = .05,
                   gralpha = .05) {
  cidf <- melt(cimat)
  names(cidf) <- c("freq", "rep", "fft")
  cigdf <- return_ci(cimat, test = "upper", alpha)
  ggplot(subset(cidf, freq <= maxfreq), aes(freq, fft, group = rep)) + geom_line(alpha =
                                                                                   gralpha) +
    geom_line(
      data = subset(cigdf, freq <= maxfreq),
      aes(freq, upper, group = 1),
      colour = "blue"
    ) +
    geom_line(
      data = subset(cigdf, freq <= maxfreq),
      aes(freq, lower, group = 1),
      colour = "blue"
    ) +
    geom_line(data = subset(estdf, freq <= maxfreq),
              colour = "red",
              aes(group = 1)) +
    theme_apa()
}

return_ci <- function(cimat,
                      test = c("lower", "two-sided", "upper"),
                      alpha = .05) {
  tst <- match.arg(test)
  ret <- data.frame(
    freq = as.numeric(rownames(cimat)),
    lower = rep(NaN, nrow(cimat)),
    upper = NaN
  )
  for (i in seq(nrow(cimat))) {
    tmp <- sort(c(cimat[i, ], -Inf, Inf))
    if (tst == "lower") {
      ret$lower[i] <- quantile(tmp, alpha)
    }
    if (tst == "upper") {
      ret$upper[i] <- quantile(tmp, 1 - alpha)
    }
    else{
      ret$lower[i] <- quantile(tmp, alpha / 2)
      ret$upper[i] <- quantile(tmp, 1 - alpha / 2)
    }
  }
  ret
}

pvals <- function(cimat,
                  estdf,
                  maxfreq = 31,
                  method = "holm") {
  ret <- data.frame(
    Frequency = as.numeric(rownames(cimat)),
    Estimate = estdf$fft,
    CIlow = -Inf,
    CIupp = Inf,
    p = NA,
    p.adj = NA
  )
  ret <- subset(ret, Frequency <= maxfreq)
  for (i in seq(nrow(ret))) {
    ## ret$CIlow[i] <- quantile(cimat[i,], .025)
    ret$CIupp[i] <- quantile(cimat[i, ], .95)
    ret$p[i] <- 1 - mean(c(cimat[i, ], Inf, -Inf) < ret$Estimate[i])
  }
  ret$p.adj <- p.adjust(ret$p, method = method)
  ret <- round(ret, 3)
  ret$p <- paste(ret$p, symnum(
    ret$p,
    corr = FALSE,
    cutpoints = c(0, .001, .01, .05, .1, 1),
    symbols = c("***", "**", "*", ".", " ")
  ))
  ret$p.adj <- paste(ret$p.adj,
                     symnum(
                       ret$p.adj,
                       corr = FALSE,
                       cutpoints = c(0, .001, .01, .05, .1, 1),
                       symbols = c("***", "**", "*", ".", " ")
                     ))
  ret
}

ci <- makedistr(df2, 10000)
plotci(ci, fftwhole(df2))

pvals(ci, fftwhole(df2))