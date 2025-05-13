## rseed <- sample(.Machine$integer.max, 1)
## set.seed(rseed)
set.seed(1601801600) # value obtained from sample(.Machine$integer.max, 1)

## devtools::install_github("bartekkroczek/detrendeR", build_vignettes=TRUE)

library(signal)
library(plyr)
library(reshape2)
library(zoo)
library(detrendeR)
library(ggplot2)
library(RColorBrewer)
library(see)

source("Symulacje/moje rekonstrukcje/functions2.r")

# Define a colorblind-friendly theme
apatheme <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12),
    legend.position = "bottom"
  )

# Option 1: Using ColorBrewer palette
cbPalette <- brewer.pal(8, "Set2")  # Set2 is colorblind friendly and prints well

## Parametry:
Np <- 1000 # liczba osób
No <- 100 # liczba obserwacji na osobę
tt <- 1 # długość interwału czasowego (s)
mnns <- .6 # minimum zakresu szumu
mxns <- .9 # maximum zakresu szumu
fr <- 6 # częstotliwość oscyl4acji (Hz) [Bartek: 6.5-7.5, 8; Chen: 5-7; Song: 3]
Tfr1 <- .5 # częstotliwość oscylacji (Hz) trendu (wolnego)
Tfr2 <- 1.5 # częstotliwość oscylacji (Hz) trendu (szybkiego)
am <- .014 # amplituda oscylacji
Tam <- .03 # amplituda trendu
padl <- 50 # długość zero-padding


# %%
# gd1 <- gen_dataset(Np, No, tt, mnns, mxns, fr, am, Tfr1, Tam)
# gd2 <- gen_dataset(Np, No, tt, mnns, mxns, fr, 0, Tfr1, Tam)
gd3 <- gen_dataset(Np, No, tt, mnns, mxns, fr, am, Tfr2, Tam)
gd4 <- gen_dataset(Np, No, tt, mnns, mxns, fr, 0, Tfr2, Tam)


# data1T <- hanning_window(preprocess(gd1, 4, padl), .(id))
# data2T <- hanning_window(preprocess(gd2, 4, padl), .(id))
data3T <- hanning_window(preprocess(gd3, 4, padl), .(id))
data4T <- hanning_window(preprocess(gd4, 4, padl), .(id))

## gd5 <- gen_dataset(Np, No, tt, mnns, mxns, fr, am*2, Tfr1, Tam)
## gd6 <- gen_dataset(Np, No, tt, mnns, mxns, fr, 0, Tfr1, Tam)
## data5T <- hanning_window(preprocess(gd5, 4, padl), .(id))
## data6T <- hanning_window(preprocess(gd6, 4, padl), .(id))

firbs <- c(0.06064464, 0.06747956, 0.07308965, 0.07725973, 0.07982848, 0.080696, 0.07982848, 0.07725973, 0.07308965, 0.06747956, 0.06064464)


dataTn <- rbind(
   # data.frame(remove_trend(data1T, type="none"),                           trtype="none",      osc="present", Tfr=Tfr1),
#   data.frame(remove_trend(data1T, type="lmf"),                            trtype="lm",        osc="present", Tfr=Tfr1),
#   data.frame(remove_trend(data1T, type="rollm", k=20),                    trtype="roll.mean", osc="present", Tfr=Tfr1),
# #  data.frame(remove_trend(data1T, type="rollr", window=No/3, step=No/27), trtype="rollreg",   osc="present", Tfr=Tfr1),
# #  data.frame(remove_trend(data1T, type="lpf", n=2, W=.1),                 trtype="lpf",       osc="present", Tfr=Tfr1),
#   data.frame(remove_trend(data1T, type="lmfnp"),                          trtype="lmf-np",    osc="present", Tfr=Tfr1),
#   data.frame(remove_trend(data1T, type="firf", n=10, w=c(0,2)/(No*tt/2)), trtype="fir",       osc="present", Tfr=Tfr1),
# #  data.frame(remove_trend(data1T, type="firf2", bs=firbs),                trtype="fir2",      osc="present", Tfr=Tfr1),
#   data.frame(remove_trend(data1T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="present", Tfr=Tfr1),
#   data.frame(remove_trend(data1T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="present", Tfr=Tfr1),
#   
#   data.frame(remove_trend(data2T, type="none"),                           trtype="none",      osc="absent", Tfr=Tfr1),
#   data.frame(remove_trend(data2T, type="lmf"),                            trtype="lm",        osc="absent", Tfr=Tfr1),
#   data.frame(remove_trend(data2T, type="rollm", k=20),                    trtype="roll.mean", osc="absent", Tfr=Tfr1),
# #  data.frame(remove_trend(data2T, type="rollr", window=No/3, step=No/27), trtype="rollreg",   osc="absent", Tfr=Tfr1),
# #  data.frame(remove_trend(data2T, type="lpf", n=2, W=.1),                 trtype="lpf",       osc="absent", Tfr=Tfr1),
#   data.frame(remove_trend(data2T, type="lmfnp"),                          trtype="lmf-np",    osc="absent", Tfr=Tfr1),
#   data.frame(remove_trend(data2T, type="firf", n=10, w=c(0,2)/(No*tt/2)), trtype="fir",       osc="absent", Tfr=Tfr1),
# #  data.frame(remove_trend(data2T, type="firf2", bs=firbs),                trtype="fir2",      osc="absent", Tfr=Tfr1),
#   data.frame(remove_trend(data2T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="absent", Tfr=Tfr1),
#   data.frame(remove_trend(data2T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="absent", Tfr=Tfr1),
  
  data.frame(remove_trend(data3T, type="none"),                           trtype="none",      osc="present", Tfr=Tfr2),
 # data.frame(remove_trend(data3T, type="lmf"),                            trtype="lm",        osc="present", Tfr=Tfr2),
###  data.frame(remove_trend(data3T, type="rollm", k=3),                    trtype="roll.mean_3", osc="present", Tfr=Tfr2),
###  data.frame(remove_trend(data3T, type="rollm", k=15),                    trtype="roll.mean_15", osc="present", Tfr=Tfr2),
###  data.frame(remove_trend(data3T, type="rollm", k=30),                    trtype="roll.mean_30", osc="present", Tfr=Tfr2),
#  data.frame(remove_trend(data3T, type="rollr", window=No/3, step=No/27), trtype="rollreg",   osc="present", Tfr=Tfr2),
#  data.frame(remove_trend(data3T, type="lpf", n=2, W=.1),                 trtype="lpf",       osc="present", Tfr=Tfr2),
 data.frame(remove_trend(data3T, type="lmfnp"),                          trtype="lmf-np",    osc="present", Tfr=Tfr2),
### data.frame(remove_trend(data3T, type="firf", n=10, w=c(0,2)/(No*tt/2)), trtype="fir",       osc="present", Tfr=Tfr2),
#  data.frame(remove_trend(data3T, type="firf2", bs=firbs),                trtype="fir2",      osc="present", Tfr=Tfr2),
#  data.frame(remove_trend(data3T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="present", Tfr=Tfr2),
 data.frame(remove_trend(data3T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="present", Tfr=Tfr2),
  
 data.frame(remove_trend(data4T, type="none"),                           trtype="none",      osc="absent", Tfr=Tfr2),
 # data.frame(remove_trend(data4T, type="lmf"),                            trtype="lm",        osc="absent", Tfr=Tfr2),
 ### data.frame(remove_trend(data4T, type="rollm", k=3),                    trtype="roll.mean_3", osc="absent", Tfr=Tfr2),
###  data.frame(remove_trend(data4T, type="rollm", k=15),                    trtype="roll.mean_15", osc="absent", Tfr=Tfr2),
###  data.frame(remove_trend(data4T, type="rollm", k=30),                    trtype="roll.mean_30", osc="absent", Tfr=Tfr2)
#  data.frame(remove_trend(data4T, type="rollr", window=No/3, step=No/27), trtype="rollreg",   osc="absent", Tfr=Tfr2),
#  data.frame(remove_trend(data4T, type="lpf", n=2, W=.1),                 trtype="lpf",       osc="absent", Tfr=Tfr2),
  data.frame(remove_trend(data4T, type="lmfnp"),                          trtype="lmf-np",    osc="absent", Tfr=Tfr2),
###  data.frame(remove_trend(data4T, type="firf", n=10, w=c(0,2)/(No*tt/2)), trtype="fir",       osc="absent", Tfr=Tfr2),
#  data.frame(remove_trend(data4T, type="firf2", bs=firbs),                trtype="fir2",      osc="absent", Tfr=Tfr2),
#  data.frame(remove_trend(data4T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="absent", Tfr=Tfr2),
  data.frame(remove_trend(data4T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="absent", Tfr=Tfr2)
)



dataTn <- rbind(
  data.frame(remove_trend(data3T, type="none"),                           trtype="none",      osc="present", Tfr=Tfr2),
  data.frame(remove_trend(data3T, type="lmfnp"),                          trtype="lmf-np",    osc="present", Tfr=Tfr2),
  data.frame(remove_trend(data3T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="present", Tfr=Tfr2),
  
  data.frame(remove_trend(data4T, type="none"),                           trtype="none",      osc="absent", Tfr=Tfr2),
  data.frame(remove_trend(data4T, type="lmfnp"),                          trtype="lmf-np",    osc="absent", Tfr=Tfr2),
  data.frame(remove_trend(data4T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="absent", Tfr=Tfr2)
)


## dataTn$samplit <- am
## dataTn <- rbind(dataTn,
##     data.frame(remove_trend(data5T, type="none"),                           trtype="none",      osc="present",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data5T, type="lmf"),                            trtype="lm",        osc="present",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data5T, type="rollm", k=20),                    trtype="roll.mean", osc="present",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data5T, type="rollr", window=No/3, step=No/27), trtype="rollreg",   osc="present",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data5T, type="lpf", n=2, W=.1),                 trtype="lpf",       osc="present",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data5T, type="lmfnp"),                          trtype="lmf-np",    osc="present",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data5T, type="firf", n=10, w=c(0,2)/(No*tt/2)), trtype="fir",       osc="present",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data5T, type="firf2", bs=firbs),                trtype="fir2",      osc="present",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data5T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="present",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data5T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="present",Tfr=Tfr1, samplit=am*2),

##     data.frame(remove_trend  apatheme
# NIe wiadomo co to robi(data6T, type="none"),                           trtype="none",      osc="absent",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data6T, type="lmf"),                            trtype="lm",        osc="absent",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data6T, type="rollm", k=20),                    trtype="roll.mean", osc="absent",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data6T, type="rollr", window=No/3, step=No/27), trtype="rollreg",   osc="absent",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data6T, type="lpf", n=2, W=.1),                 trtype="lpf",       osc="absent",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data6T, type="lmfnp"),                          trtype="lmf-np",    osc="absent",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data6T, type="firf", n=10, w=c(0,2)/(No*tt/2)), trtype="fir",       osc="absent",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data6T, type="firf2", bs=firbs),                trtype="fir2",      osc="absent",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data6T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="absent",Tfr=Tfr1, samplit=am*2),
##     data.frame(remove_trend(data6T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="absent",Tfr=Tfr1, samplit=am*2)
## )

dataTn$len <- tt
dataTn <- remove_padding(dataTn)
#dataTn$trtype <- factor(dataTn$trtype,
#                         levels = c("fir", "lmf-np", "none", "polyn-np", "roll.mean"),
#                         labels = c("FIR", "LM", "Raw", "Polynomial", "Rolling Mean"))

dataTn$trtype <- factor(dataTn$trtype,
                         levels = c("none", "lmf-np", "polyn-np"),
                         labels = c("Raw", "Linear", "Polynomial deg = 4"))
ggplot(aggregate(trend ~ time * trtype * osc * Tfr, dataTn, FUN = mean), 
       aes(time, trend)) +
  geom_line(size = 1.1, 
            alpha = 0.7) +
  facet_grid(interaction(osc, Tfr) ~ trtype,
             labeller = label_wrap_gen(width = 10)) +  # Wrap long facet labels
  scale_x_continuous("Time") +
  scale_y_continuous("Trend") +
  apatheme

dataFn <- fourier(dataTn, .(id, trtype, osc, Tfr))
dataFna <- subset(aggregate(fft ~ freq * trtype * osc * Tfr, dataFn, FUN = mean), freq <= No / 2)
# dataFna$trtype <- factor(dataFna$trtype,
#                         levels = c("fir", "lmf-np", "none", "polyn-np", "roll.mean"),
#                         labels = c("FIR", "LM", "Raw", "Polynomial", "Rolling Mean"))


# dataFn <- fourier(dataTn, .(id, trtype, osc, Tfr, samplit))
# dataFna <- subset(aggregate(fft~freq*trtype*osc*samplit, subset(dataFn, Tfr==Tfr1), FUN=mean), freq<=No/2)
# ggplot(subset(dataFna, freq<=32), aes(freq, fft, colour=trtype, linetype=trtype)) + geom_line() + facet_grid(factor(samplit)~osc)
# %% TO SĄ WYKRESY SYMULACJI!!!!!!

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
okabe_ito <- c("#E69F00", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#000000")
# subset(subset(dataFna, freq <= 20), trtype %in% c("Raw", "FIR", "Rolling Mean"))
ggplot(subset(dataFna, freq <= 20), aes(freq, fft, colour = trtype)) +
  geom_line(size = 1.1, alpha = 0.7) +  # Added alpha for transparency
  facet_wrap(~ osc) +  # Removed Tfr from faceting
  scale_x_continuous("Frequency") +
  scale_y_continuous("Power") +
  scale_colour_manual("Detrending method", values = okabe_ito) +  # Changed scale_fill_brewer to scale_fill_manual
  geom_text(
    data = data.frame(
      x = c(1, 4.8,  9.2, 1, 6, 10),
      y = c(1.2, .92, .92,  1.2, 1.0, 0.95),
      label = c("A", "", "", "A", "B", ""),
      osc = c("absent", "absent","absent",  "present", "present", "present")
    ),
    aes(x, y, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  apatheme
# NIe wiadomo co to robi

# 
# ci1 <- makedistr(gd2, 1000, FUN = function(data) {
#   aggregate(fft ~ freq,
#     FUN = mean,
#     fourier(
#       remove_padding(data.frame(
#         remove_trend(
#           hanning_window(
#             preprocess(data, 4, padl), .(id)
#           ),
#           type = "firf", n = 10, w = c(0, 2) / (No * tt / 2)
#         )
#       )), .(id)
#     )
#   )
# })
# # %%
# 
# plotci(ci1, aggregate(fft ~ freq, FUN = mean, subset(dataFn, trtype == "fir" & osc == "absent" & Tfr == Tfr1)))
# pvals(ci1, aggregate(fft ~ freq, FUN = mean, subset(dataFn, trtype == "fir" & osc == "absent" & Tfr == Tfr1)))
# 
# # %%
# 
# ci2 <- makedistr(gd1, 1000, FUN = function(data) {
#   aggregate(fft ~ freq,
#     FUN = mean,
#     fourier(
#       remove_padding(data.frame(
#         remove_trend(
#           hanning_window(
#             preprocess(data, 4, padl), .(id)
#           ),
#           type = "firf", n = 10, w = c(0, 2) / (No * tt / 2)
#         )
#       )), .(id)
#     )
#   )
# })
## plotci(ci2, aggregate(fft~freq, FUN=mean, subset(dataFn, trtype=="fir" & osc=="present" & Tfr==Tfr1)))
## pvals(ci2, aggregate(fft~freq, FUN=mean, subset(dataFn, trtype=="fir" & osc=="present" & Tfr==Tfr1)))
# %%
## ci3 <- makedistr(gd2, 1000, FUN = function(data){
##   aggregate(fft~freq, FUN=mean,
##             fourier(
##                 remove_padding(data.frame(
##                     remove_trend(
##                         hanning_window(
##                             preprocess(data, 4, padl), .(id)), type="polnp", degree=2*2*tt))), .(id)))})
## plotci(ci3, aggregate(fft~freq, FUN=mean, subset(dataFn, trtype=="polyn-np" & osc=="absent" & Tfr==Tfr1)))
## pvals(ci3, aggregate(fft~freq, FUN=mean, subset(dataFn, trtype=="polyn-np" & osc=="absent" & Tfr==Tfr1)))

## ci4 <- makedistr(gd1, 1000, FUN = function(data){
##   aggregate(fft~freq, FUN=mean,
##             fourier(
##                 remove_padding(data.frame(
##                     remove_trend(
##                         hanning_window(
##                             preprocess(data, 4, padl), .(id)), type="polnp", degree=2*2*tt))), .(id)))})
## plotci(ci4, aggregate(fft~freq, FUN=mean, subset(dataFn, trtype=="polyn-np" & osc=="absent" & Tfr==Tfr1)))
## pvals(ci4, aggregate(fft~freq, FUN=mean, subset(dataFn, trtype=="polyn-np" & osc=="absent" & Tfr==Tfr1)))

makedistr







# %%

plotci <- function(cimat,
                   estdf,
                   maxfreq = 31,
                   alpha = .05,
                   gralpha = .05) {
  cidf <- melt(cimat)
  names(cidf) <- c("freq", "rep", "fft")
  cigdf <- return_ci(cimat, test = "upper", alpha)
  ggplot(subset(cidf, freq <= maxfreq), aes(freq, fft, group = rep)) +
    geom_line(
      alpha =
        gralpha
    ) +
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
    geom_line(
      data = subset(estdf, freq <= maxfreq),
      colour = "red",
      aes(group = 1)
    )
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
    } else {
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
  ret$p.adj <- paste(
    ret$p.adj,
    symnum(
      ret$p.adj,
      corr = FALSE,
      cutpoints = c(0, .001, .01, .05, .1, 1),
      symbols = c("***", "**", "*", ".", " ")
    )
  )
  ret
}


fftwhole <- function(data, FUN = mean) {
  temp <- remove_trend(
    # preprocess removing only very extreame values and no padding at all
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
    {
      (. - 1) / diff(range(data$time))
    }

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


# %%
