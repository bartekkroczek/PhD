#%% 
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

print(getwd())
source("Symulacje/moje rekonstrukcje/functions2.r")
#%% 
#repo_low <- "https://raw.githubusercontent.com/bartekkroczek/Behavioral-oscillations/master/data/low_resolution/"
#files_low <- read.table("files_low.txt")$V1
freq_low <- 1 / 60

# d1 <- NULL
#for (i in files_low)
#  d1 <- rbind(d1, read.csv(paste0(repo_low, i)))
#d1$Corr = as.numeric(as.factor(d1$Corr)) - 1

csv_files <- dir_ls(path = "./Experiment_1/Data/", glob = "*.csv")

d1 <- csv_files %>%
  map_df(~read_csv(.))


d1 <- dplyr::rename(d1, id = "PART_ID")
d1 <- dplyr::rename(d1, corr = "Corr")
d1 <- subset(d1, Trial_type == "experiment")
d1atemp <- aggregate(corr ~ id, d1, FUN = mean)
d1a <- aggregate(corr ~ CSI * id, d1, FUN = mean)
d1a <- subset(d1a, id %in% d1atemp$id[d1atemp$corr > .6])
d1a$time <- with(d1a, CSI * freq_low)

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
              aes(group = 1))
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

ci <- makedistr(d1a, 10)
plotci(ci, fftwhole(d1a))

pvals(ci, fftwhole(d1a))