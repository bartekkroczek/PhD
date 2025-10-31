# install.packages("multitaper")  # uncomment if needed
library(multitaper)

set.seed(42)

# 1) Toy data: three short series
N  <- 256       # sample count (short)
fs <- 1         # Hz (sampling rate)
dt <- 1/fs
t  <- seq(0, (N-1)/fs, by=dt)

# Series A: clear 0.15 Hz sinusoid + noise
fA <- 0.15
xA <- 1.5*sin(2*pi*fA*t + 0.5) + rnorm(N, 0, 1.0)

# Series B: two sinusoids (0.28 Hz, 0.40 Hz) + noise
fB1 <- 0.28; fB2 <- 0.40
xB <- 1.0*sin(2*pi*fB1*t) + 0.8*sin(2*pi*fB2*t + 1.2) + rnorm(N, 0, 1.2)

# Series C: pure noise (null)
xC <- rnorm(N, 0, 1.0)

# 2) MTM parameters: choose NW and K
NW <- 3
K  <- 2*NW - 1  # common choice; you can also use 2*NW
# -> spectral half-bandwidth W = NW/T, where T = N*dt

# 3) Compute MTM spectra + harmonic F-test
mtmA <- spec.mtm(ts(xA, deltat=dt), nw=NW, k=K, Ftest=TRUE, plot=FALSE, deltat=dt, dtUnits="s")
mtmB <- spec.mtm(ts(xB, deltat=dt), nw=NW, k=K, Ftest=TRUE, plot=FALSE, deltat=dt, dtUnits="s")
mtmC <- spec.mtm(ts(xC, deltat=dt), nw=NW, k=K, Ftest=TRUE, plot=FALSE, deltat=dt, dtUnits="s")

# Helper to safely extract F-test vector from 'mtm' objects across versions
get_F <- function(obj) {
  # Try common slots used by multitaper::spec.mtm over versions
  if (!is.null(obj$mtm$Ftest)) return(obj$mtm$Ftest)
  if (!is.null(obj$Ftest))     return(obj$Ftest)
  if (!is.null(obj$mtm$mtm$Ftest)) return(obj$mtm$mtm$Ftest)
  stop("Could not locate F-test values in the mtm object. Run str(obj) to inspect slots.")
}

fAvals <- get_F(mtmA)
fBvals <- get_F(mtmB)
fCvals <- get_F(mtmC)

freqA <- mtmA$freq
freqB <- mtmB$freq
freqC <- mtmC$freq

# 4) Visualization: spectrum (top row) and F-test (bottom row) with significance lines
#    Under H0, F ~ F(2, 2K-2). Choose e.g. 0.95 and 0.99 significance.
sig_lines <- c(0.95, 0.99)
df1 <- 2
df2 <- 2*K - 2
thresh <- qf(sig_lines, df1, df2)

par(mfrow=c(2,3), mar=c(4,4,2,1))

# Top row: spectra (using built-in plot.mtm for consistency)
plot(mtmA, Ftest=FALSE, main="Series A: MTM Spectrum")
abline(v=fA, col="forestgreen", lty=2)

plot(mtmB, Ftest=FALSE, main="Series B: MTM Spectrum")
abline(v=c(fB1, fB2), col="forestgreen", lty=2)

plot(mtmC, Ftest=FALSE, main="Series C: MTM Spectrum")

# Bottom row: harmonic F-test with significance lines
plot(mtmA, Ftest=TRUE, main="Series A: Harmonic F-test")
abline(h=thresh, col=c("orange","red"), lty=3)
abline(v=fA, col="forestgreen", lty=2)
legend("topright", legend=paste0("p=", sig_lines), lty=3, col=c("orange","red"), bty="n")

plot(mtmB, Ftest=TRUE, main="Series B: Harmonic F-test")
abline(h=thresh, col=c("orange","red"), lty=3)
abline(v=c(fB1, fB2), col="forestgreen", lty=2)

plot(mtmC, Ftest=TRUE, main="Series C: Harmonic F-test")
abline(h=thresh, col=c("orange","red"), lty=3)

# 5) Mark statistically significant local peaks (e.g., p <= 0.01) on F-test plots
#    A simple local-maximum finder:
local_maxima <- function(y) {
  # returns indices of strict local maxima
  dy <- diff(y)
  which(diff(sign(dy)) == -2) + 1
}

p_cut <- 0.01
F_cut <- qf(1 - p_cut, df1, df2)

mark_peaks <- function(freq, Fvals, label="", col="red") {
  idx <- local_maxima(Fvals)
  idx_sig <- idx[Fvals[idx] >= F_cut]
  points(freq[idx_sig], Fvals[idx_sig], pch=19, col=col)
  if (length(idx_sig) > 0) {
    mtext(paste0(label, " sig peaks: ",
                 paste(round(freq[idx_sig], 3), collapse=", ")),
          side=3, line=-1.2, adj=0, cex=0.7, col=col)
  }
}

# Re-plot the bottom row for marking (or run the plotting above again)
par(mfrow=c(1,3), mar=c(4,4,2,1))
plot(mtmA, Ftest=TRUE, main="A: F-test (sig peaks in red)")
abline(h=F_cut, col="red", lty=3)
mark_peaks(freqA, fAvals, "A")
abline(v=fA, col="forestgreen", lty=2)

plot(mtmB, Ftest=TRUE, main="B: F-test (sig peaks in red)")
abline(h=F_cut, col="red", lty=3)
mark_peaks(freqB, fBvals, "B")
abline(v=c(fB1, fB2), col="forestgreen", lty=2)

plot(mtmC, Ftest=TRUE, main="C: F-test (sig peaks in red)")
abline(h=F_cut, col="red", lty=3)
mark_peaks(freqC, fCvals, "C")

