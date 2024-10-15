library(imputeTS)

gen_signal <- function(Npart, Nobs, time, freq, ampl) {
  ret <- NULL
  for (i in seq(Npart)) {
    phase <- phase <- runif(1, 0, 2 * pi)
    ret <- c(ret, sin(seq(Nobs) / Nobs * (2 * pi * time * freq) + phase) *
               ampl)
  }
  ret
}

gen_dataset <- function(Npart,
                        Nobs,
                        time,
                        mnns,
                        mxns,
                        freq,
                        ampl,
                        trendfreq,
                        trendampl) {
  data.frame(
    id = paste0("p", rep(seq(Npart), each = Nobs)),
    time = rep(seq(0, time, length.out = Nobs), Npart),
    corr = runif(Npart * Nobs, mnns, mxns) +
      gen_signal(Npart, Nobs, time, freq, ampl) +
      gen_signal(Npart, Nobs, time, trendfreq, trendampl)
  )
}

preprocess <- function(data,
                       outlier = Inf,
                       padlen = 50) {
  ## normalization
  ret <- ddply(data, .(id), function(df) {
    cbind(df$time, melt(as.vector(scale(df$corr))))
  })
  ret <- dplyr::rename(ret, corr = "value", time = "df$time")
  ## outlier removal
  ret <- subset(ret, abs(corr) < outlier)
  ## zero padding
  ret$padding <- 0
  timestep <- diff(ret$time)
  timestep <- mean(timestep[timestep > 0])
  mit <- min(ret$time)
  mat <- max(ret$time)
  if (padlen > 0)
    for (i in unique(ret$id))
      ret <- rbind(ret, data.frame(
        id = i,
        time = c(
          seq(
            to = mit - timestep,
            length.out = padlen,
            by = timestep
          ),
          seq(
            from = mat + timestep,
            length.out = padlen,
            by = timestep
          )
        ),
        corr = 0,
        padding = 1
      ))
  ret <- ret[order(ret$time), ]
  ret[order(ret$id), ]
}

hanning_window <- function(data, columns) {
  ret <- ddply(data, columns, function(df) {
    cbind(df$time, melt(df$corr * hanning(length(df$corr))), df$padding)
  })
  dplyr::rename(ret,
                ,
                time = "df$time",
                corr = "value",
                padding = "df$padding")
}

remove_trend <- function(data,
                         type = c(
                           "none",
                           "meanf",
                           "lmf",
                           "rollm",
                           "rollr",
                           "lpf",
                           "lmfnp",
                           "firf",
                           "pol",
                           "firf2",
                           "polnp"
                         ),
                         ...) {
  type <- match.arg(type)
  
  none <- function(dtf, ...)
    melt(rep(0, nrow(dtf)))
  meanf <- function(dtf, ...)
    melt(rep(mean(dtf$corr), nrow(dtf)))
  lmf <- function(dtf, ...)
    melt(as.vector(getltrend(dtf$corr, dtf$time)))
  lmfnp <- function(dtf, ...)
    melt(as.vector(predict(
      with(subset(dtf, padding == 0), getltrend(corr, time)), newtime = dtf$time
    )))
  rollm <- function(dtf, ...)
    melt(rollmean(dtf$corr, na.pad = TRUE, ...))
  rollr <- function(dtf, ...)
    melt(rollregres(dtf$corr, ...))
  lpf <- function(dtf, ...)
    melt(as.vector(signal::filter(butter(type = "low", ...), dtf$corr)))
  firf <- function(dtf, ...)
    melt(as.vector(signal::filter(
      fir1(type = "pass", window = boxcar, ...), dtf$corr
    )))
  firf2 <- function(dtf, ...)
    melt(FIR(dtf$corr, ...))
  pol <- function(dtf, ...)
    melt(as.vector(gettrend(dtf$corr, time = dtf$time, ...)))
  polnp <- function(dtf, ...)
    melt(as.vector(predict(
      with(subset(dtf, padding == 0), gettrend(corr, time = time, ...)), newtime =
        dtf$time
    )))
  
  data$trend <- ddply(data, .(id), get(type), ...)$value
  
  data$corr <- data$corr - data$trend
  data
}

remove_padding <- function(data) {
  subset(data, padding == 0)
}

FIR <- function(signal, bs) {
  lb <- length(bs)
  if (lb > length(signal))
    stop("The signal is to short for the given impulse response vector.")
  ret <- rep(NA, lb)
  for (i in lb:length(signal)) {
    ret[i] <- sum(bs * signal[(i - lb + 1):i])
  }
  ret
}

rollregres <- function(signal, window, step) {
  window <- round(window)
  step <- round(step)
  intcp <- NULL
  slope <- NULL
  for (i in seq(1, length(signal), step)) {
    if (i < window) {
      intcp[i] <- NA
      slope[i] <- NA
    }
    else{
      windex <- seq(i - window + 1, i)
      m <- lm(signal[windex] ~ windex)
      intcp[i] <- coef(m)[1]
      slope[i] <- coef(m)[2]
    }
  }
  coefs <- data.frame(
    intercept = na_ma(intcp, k = 1, weighting = "simple"),
    slope = na_ma(slope, k = 1, weighting = "simple")
  )
  pred <- NULL
  for (i in seq(length(signal)))
    pred[i] <- coefs$intercept[i] + coefs$slope[i] * i
  pred
}

fourier <- function(data, columns) {
  data <- subset(data, !is.na(data$corr))
  p1 <- ddply(data, columns, function(df)
    melt(Mod(fft(df$corr) / sqrt(nrow(
      df
    )))))
  p1 <- dplyr::rename(p1, fft = "value")
  p2 <- ddply(data, columns, function(df)
    melt((seq(nrow(
      df
    )) - 1) / df$len))
  data.frame(p1, freq = p2$value)
}
