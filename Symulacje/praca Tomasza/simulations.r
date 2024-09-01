set.seed(1601801600) # value obtained from sample(.Machine$integer.max, 1)

library(zoo)

source("functions.r")

## Parametry:
Np   <- 1000  # number of agents
No   <- 100   # samples per agent
tt   <- 1     # time inteval (s)
mnns <- .6    # minimal noise
mxns <- .9    # maximal noise
fr   <- 6     # oscilations frequency (Hz)
Tfr1 <- .5    # (slow) trend frequency (Hz)
Tfr2 <- 1.5   # (fast) trend frequency (Hz)
am   <- .014  # oscilations amplitude
Tam  <- .03   # trend amplitude
padl <- 50    # zer-padding lenght (in samples)

gd1 <- gen_dataset(Np, No, tt, mnns, mxns, fr, am, Tfr1, Tam)
gd2 <- gen_dataset(Np, No, tt, mnns, mxns, fr,  0, Tfr1, Tam)
gd3 <- gen_dataset(Np, No, tt, mnns, mxns, fr, am, Tfr2, Tam)
gd4 <- gen_dataset(Np, No, tt, mnns, mxns, fr,  0, Tfr2, Tam)

data1T <- hanning_window(preprocess(gd1, 4, padl), .(id))
data2T <- hanning_window(preprocess(gd2, 4, padl), .(id))
data3T <- hanning_window(preprocess(gd3, 4, padl), .(id))
data4T <- hanning_window(preprocess(gd4, 4, padl), .(id))

#### dataTn <- rbind(
####     data.frame(remove_trend(data1T, type="none"),                           trtype="none",      osc="present", Tfr=Tfr1),
####     data.frame(remove_trend(data1T, type="lmf"),                            trtype="lm",        osc="present", Tfr=Tfr1),
####     data.frame(remove_trend(data1T, type="rollm", k=20),                    trtype="roll.mean", osc="present", Tfr=Tfr1),
####     data.frame(remove_trend(data1T, type="lmfnp"),                          trtype="lmf-np",    osc="present", Tfr=Tfr1),
####     data.frame(remove_trend(data1T, type="firf2", bs=firbs),                trtype="fir2",      osc="present", Tfr=Tfr1),
####     data.frame(remove_trend(data1T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="present", Tfr=Tfr1),
####     data.frame(remove_trend(data1T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="present", Tfr=Tfr1),

####     data.frame(remove_trend(data2T, type="none"),                           trtype="none",      osc="absent", Tfr=Tfr1),
####     data.frame(remove_trend(data2T, type="lmf"),                            trtype="lm",        osc="absent", Tfr=Tfr1),
####     data.frame(remove_trend(data2T, type="rollm", k=20),                    trtype="roll.mean", osc="absent", Tfr=Tfr1),
####     data.frame(remove_trend(data2T, type="lmfnp"),                          trtype="lmf-np",    osc="absent", Tfr=Tfr1),
####     data.frame(remove_trend(data2T, type="firf2", bs=firbs),                trtype="fir2",      osc="absent", Tfr=Tfr1),
####     data.frame(remove_trend(data2T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="absent", Tfr=Tfr1),
####     data.frame(remove_trend(data2T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="absent", Tfr=Tfr1),

####     data.frame(remove_trend(data3T, type="none"),                           trtype="none",      osc="present", Tfr=Tfr2),
####     data.frame(remove_trend(data3T, type="lmf"),                            trtype="lm",        osc="present", Tfr=Tfr2),
####     data.frame(remove_trend(data3T, type="rollm", k=20),                    trtype="roll.mean", osc="present", Tfr=Tfr2),
####     data.frame(remove_trend(data3T, type="lmfnp"),                          trtype="lmf-np",    osc="present", Tfr=Tfr2),
####     data.frame(remove_trend(data3T, type="firf2", bs=firbs),                trtype="fir2",      osc="present", Tfr=Tfr2),
####     data.frame(remove_trend(data3T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="present", Tfr=Tfr2),
####     data.frame(remove_trend(data3T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="present", Tfr=Tfr2),

####     data.frame(remove_trend(data4T, type="none"),                           trtype="none",      osc="absent", Tfr=Tfr2),
####     data.frame(remove_trend(data4T, type="lmf"),                            trtype="lm",        osc="absent", Tfr=Tfr2),
####     data.frame(remove_trend(data4T, type="rollm", k=20),                    trtype="roll.mean", osc="absent", Tfr=Tfr2),
####     data.frame(remove_trend(data4T, type="lmfnp"),                          trtype="lmf-np",    osc="absent", Tfr=Tfr2),
####     data.frame(remove_trend(data4T, type="firf2", bs=firbs),                trtype="fir2",      osc="absent", Tfr=Tfr2),
####     data.frame(remove_trend(data4T, type="pol", degree=2*2*tt),             trtype="polyn",     osc="absent", Tfr=Tfr2),
####     data.frame(remove_trend(data4T, type="polnp", degree=2*2*tt),           trtype="polyn-np",  osc="absent", Tfr=Tfr2)
#### )

dataTn$len=tt
dataTn <- remove_padding(dataTn)
dataFn <- fourier(dataTn, .(id, trtype, osc, Tfr))
dataFna <- subset(aggregate(fft~freq*trtype*osc*Tfr, dataFn, FUN=mean), freq<=No/2)
