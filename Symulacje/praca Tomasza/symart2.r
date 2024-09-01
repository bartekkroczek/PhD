## rseed <- sample(.Machine$integer.max, 1)
## set.seed(rseed)
set.seed(1601801600) # value obtained from sample(.Machine$integer.max, 1)

## devtools::install_github("bartekkroczek/detrendeR", build_vignettes=TRUE)

library(signal)
library(plyr)
library(reshape2)
library(zoo)
library(detrendeR)

source("functions2.r")

## Parametry:
Np   <- 1000  # liczba osób
No   <- 100   # liczba obserwacji na osobę
tt   <- 1     # długość interwału czasowego (s)
mnns <- .6    # minimum zakresu szumu
mxns <- .9    # maximum zakresu szumu
fr   <- 6     # częstotliwość oscylacji (Hz) [Bartek: 6.5-7.5, 8; Chen: 5-7; Song: 3]
Tfr1 <- .5    # częstotliwość oscylacji (Hz) trendu (wolnego)
Tfr2 <- 1.5   # częstotliwość oscylacji (Hz) trendu (szybkiego)
am   <- .014  # amplituda oscylacji
Tam  <- .03   # amplituda trendu
padl <- 50    # długość zero-padding

gd1 <- gen_dataset(Np, No, tt, mnns, mxns, fr, am, Tfr1, Tam)
gd2 <- gen_dataset(Np, No, tt, mnns, mxns, fr,  0, Tfr1, Tam)
gd3 <- gen_dataset(Np, No, tt, mnns, mxns, fr, am, Tfr2, Tam)
gd4 <- gen_dataset(Np, No, tt, mnns, mxns, fr,  0, Tfr2, Tam)


data1T <- hanning_window(preprocess(gd1, 4, padl), .(id))
data2T <- hanning_window(preprocess(gd2, 4, padl), .(id))
data3T <- hanning_window(preprocess(gd3, 4, padl), .(id))
data4T <- hanning_window(preprocess(gd4, 4, padl), .(id))

## gd5 <- gen_dataset(Np, No, tt, mnns, mxns, fr, am*2, Tfr1, Tam)
## gd6 <- gen_dataset(Np, No, tt, mnns, mxns, fr, 0, Tfr1, Tam)
## data5T <- hanning_window(preprocess(gd5, 4, padl), .(id))
## data6T <- hanning_window(preprocess(gd6, 4, padl), .(id))

firbs <- c(0.06064464, 0.06747956, 0.07308965, 0.07725973, 0.07982848, 0.080696, 0.07982848, 0.07725973, 0.07308965, 0.06747956, 0.06064464)

dataTn <- rbind(
    data.frame(remove_trend(data1T, type="none"),                           trtype="none",      osc="present", Tfr=Tfr1),
    data.frame(remove_trend(data1T, type="lmf"),                            trtype="lm",        osc="present", Tfr=Tfr1),
    data.frame(remove_trend(data1T, type="rollm", k=20),                    trtype="roll.mean", osc="present", Tfr=Tfr1),
    data.frame(remove_trend(data1T, type="rollr", window=No/3, step=No/27), trtype="rollreg",   osc="present", Tfr=Tfr1),
    data.frame(remove_trend(data1T, type="lpf", n=2, W=.1),                 trtype="lpf",       osc="present", Tfr=Tfr1),
    data.frame(remove_trend(data1T, type="lmfnp"),                          trtype="lmf-np",    osc="present", Tfr=Tfr1),
    data.frame(remove_trend(data1T, type="firf", n=10, w=c(0,2)/(No*tt/2)), trtype="fir",       osc="present", Tfr=Tfr1),
    data.frame(remove_trend(data1T, type="firf2", bs=firbs),                trtype="fir2",      osc="present", Tfr=Tfr1),
    data.frame(remove_trend(data1T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="present", Tfr=Tfr1),
    data.frame(remove_trend(data1T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="present", Tfr=Tfr1),

    data.frame(remove_trend(data2T, type="none"),                           trtype="none",      osc="absent", Tfr=Tfr1),
    data.frame(remove_trend(data2T, type="lmf"),                            trtype="lm",        osc="absent", Tfr=Tfr1),
    data.frame(remove_trend(data2T, type="rollm", k=20),                    trtype="roll.mean", osc="absent", Tfr=Tfr1),
    data.frame(remove_trend(data2T, type="rollr", window=No/3, step=No/27), trtype="rollreg",   osc="absent", Tfr=Tfr1),
    data.frame(remove_trend(data2T, type="lpf", n=2, W=.1),                 trtype="lpf",       osc="absent", Tfr=Tfr1),
    data.frame(remove_trend(data2T, type="lmfnp"),                          trtype="lmf-np",    osc="absent", Tfr=Tfr1),
    data.frame(remove_trend(data2T, type="firf", n=10, w=c(0,2)/(No*tt/2)), trtype="fir",       osc="absent", Tfr=Tfr1),
    data.frame(remove_trend(data2T, type="firf2", bs=firbs),                trtype="fir2",      osc="absent", Tfr=Tfr1),
    data.frame(remove_trend(data2T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="absent", Tfr=Tfr1),
    data.frame(remove_trend(data2T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="absent", Tfr=Tfr1),

    data.frame(remove_trend(data3T, type="none"),                           trtype="none",      osc="present", Tfr=Tfr2),
    data.frame(remove_trend(data3T, type="lmf"),                            trtype="lm",        osc="present", Tfr=Tfr2),
    data.frame(remove_trend(data3T, type="rollm", k=20),                    trtype="roll.mean", osc="present", Tfr=Tfr2),
    data.frame(remove_trend(data3T, type="rollr", window=No/3, step=No/27), trtype="rollreg",   osc="present", Tfr=Tfr2),
    data.frame(remove_trend(data3T, type="lpf", n=2, W=.1),                 trtype="lpf",       osc="present", Tfr=Tfr2),
    data.frame(remove_trend(data3T, type="lmfnp"),                          trtype="lmf-np",    osc="present", Tfr=Tfr2),
    data.frame(remove_trend(data3T, type="firf", n=10, w=c(0,2)/(No*tt/2)), trtype="fir",       osc="present", Tfr=Tfr2),
    data.frame(remove_trend(data3T, type="firf2", bs=firbs),                trtype="fir2",      osc="present", Tfr=Tfr2),
    data.frame(remove_trend(data3T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="present", Tfr=Tfr2),
    data.frame(remove_trend(data3T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="present", Tfr=Tfr2),

    data.frame(remove_trend(data4T, type="none"),                           trtype="none",      osc="absent", Tfr=Tfr2),
    data.frame(remove_trend(data4T, type="lmf"),                            trtype="lm",        osc="absent", Tfr=Tfr2),
    data.frame(remove_trend(data4T, type="rollm", k=20),                    trtype="roll.mean", osc="absent", Tfr=Tfr2),
    data.frame(remove_trend(data4T, type="rollr", window=No/3, step=No/27), trtype="rollreg",   osc="absent", Tfr=Tfr2),
    data.frame(remove_trend(data4T, type="lpf", n=2, W=.1),                 trtype="lpf",       osc="absent", Tfr=Tfr2),
    data.frame(remove_trend(data4T, type="lmfnp"),                          trtype="lmf-np",    osc="absent", Tfr=Tfr2),
    data.frame(remove_trend(data4T, type="firf", n=10, w=c(0,2)/(No*tt/2)), trtype="fir",       osc="absent", Tfr=Tfr2),
    data.frame(remove_trend(data4T, type="firf2", bs=firbs),                trtype="fir2",      osc="absent", Tfr=Tfr2),
    data.frame(remove_trend(data4T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="absent", Tfr=Tfr2),
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

##     data.frame(remove_trend(data6T, type="none"),                           trtype="none",      osc="absent",Tfr=Tfr1, samplit=am*2),
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

dataTn$len=tt
dataTn <- remove_padding(dataTn)

#### ggplot(aggregate(trend~time*trtype*osc*Tfr, dataTn, FUN=mean), aes(time, trend)) + geom_line() + facet_grid(interaction(osc,Tfr)~trtype)

dataFn <- fourier(dataTn, .(id, trtype, osc, Tfr))
dataFna <- subset(aggregate(fft~freq*trtype*osc*Tfr, dataFn, FUN=mean), freq<=No/2)

## dataFn <- fourier(dataTn, .(id, trtype, osc, Tfr, samplit))
## dataFna <- subset(aggregate(fft~freq*trtype*osc*samplit, subset(dataFn, Tfr==Tfr1), FUN=mean), freq<=No/2)
## ggplot(subset(dataFna, freq<=32), aes(freq, fft, colour=trtype, linetype=trtype)) + geom_line() + facet_grid(factor(samplit)~osc)

#### ggplot(subset(dataFna, freq<=32), aes(freq, fft, colour=trtype, linetype=trtype)) + geom_line() + facet_grid(Tfr~osc)
#### ggplot(subset(dataFna, freq<=32), aes(freq, fft, colour=trtype)) + geom_line() + facet_wrap(Tfr~osc) + scale_x_continuous("Częstotliwość") + scale_y_continuous("Moc") + scale_colour_discrete("Trend") + geom_text(data=data.frame(x=c(1,7,1,6), y=c(1.3,1,1.3,1.55), label=c("A","B","A","C"), osc=c("absent","absent","present","present")), aes(x,y,label=label), inherit.aes=FALSE, size=6)

## ci1 <- makedistr(gd2, 1000, FUN = function(data){
##   aggregate(fft~freq, FUN=mean,
##             fourier(
##                 remove_padding(data.frame(
##                     remove_trend(
##                         hanning_window(
##                             preprocess(data, 4, padl), .(id)), type="firf", n=10, w=c(0,2)/(No*tt/2)))), .(id)))})
## plotci(ci1, aggregate(fft~freq, FUN=mean, subset(dataFn, trtype=="fir" & osc=="absent" & Tfr==Tfr1)))
## pvals(ci1, aggregate(fft~freq, FUN=mean, subset(dataFn, trtype=="fir" & osc=="absent" & Tfr==Tfr1)))

## ci2 <- makedistr(gd1, 1000, FUN = function(data){
##   aggregate(fft~freq, FUN=mean,
##             fourier(
##                 remove_padding(data.frame(
##                     remove_trend(
##                         hanning_window(
##                             preprocess(data, 4, padl), .(id)), type="firf", n=10, w=c(0,2)/(No*tt/2)))), .(id)))})
## plotci(ci2, aggregate(fft~freq, FUN=mean, subset(dataFn, trtype=="fir" & osc=="present" & Tfr==Tfr1)))
## pvals(ci2, aggregate(fft~freq, FUN=mean, subset(dataFn, trtype=="fir" & osc=="present" & Tfr==Tfr1)))

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
