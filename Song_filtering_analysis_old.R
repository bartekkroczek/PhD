library(tidyverse)
library(signal)
library(pracma)
library(patchwork)

# EEGLAB-style FIR filter design
design_eeglab_fir <- function(low_freq, high_freq, fs, order = 10) {
  # Following EEGLAB's pop_eegfilt() legacy implementation
  
  nyquist <- fs/2
  if (is.null(low_freq)) low_freq <- 0
  if (is.null(high_freq)) high_freq <- nyquist
  
  # EEGLAB uses different transition bandwidths based on frequency
  if (low_freq < 1) {
    trans_band_low <- 0.5   # 0.5 Hz transition for very low frequencies
  } else if (low_freq < 2) {
    trans_band_low <- 0.75
  } else if (low_freq < 5) {
    trans_band_low <- 1
  } else if (low_freq < 10) {
    trans_band_low <- 1.5
  } else {
    trans_band_low <- 2
  }
  
  if (high_freq < 1) {
    trans_band_high <- 0.5
  } else if (high_freq < 2) {
    trans_band_high <- 0.75
  } else if (high_freq < 5) {
    trans_band_high <- 1
  } else if (high_freq < 10) {
    trans_band_high <- 1.5
  } else {
    trans_band_high <- 2
  }
  
  # Calculate filter frequencies and amplitudes
  if (low_freq <= 0) {  # Lowpass only
    freq_vector <- c(0, 
                     high_freq - trans_band_high, 
                     high_freq + trans_band_high, 
                     nyquist)
    amp_vector <- c(1, 1, 0, 0)
  } else if (high_freq >= nyquist) {  # Highpass only
    freq_vector <- c(0,
                     low_freq - trans_band_low,
                     low_freq + trans_band_low,
                     nyquist)
    amp_vector <- c(0, 0, 1, 1)
  } else {  # Bandpass
    freq_vector <- c(0,
                     low_freq - trans_band_low,
                     low_freq + trans_band_low,
                     high_freq - trans_band_high,
                     high_freq + trans_band_high,
                     nyquist)
    amp_vector <- c(0, 0, 1, 1, 0, 0)
  }
  
  # Normalize frequencies
  freq_vector <- freq_vector / nyquist
  
  # Design filter using window method (EEGLAB legacy approach)
  n_fft <- order + 1
  
  # Create ideal frequency response
  f <- seq(0, 1, length.out = n_fft)
  H <- numeric(n_fft)
  
  # Interpolate amplitude response
  for (i in 1:(length(freq_vector)-1)) {
    idx <- f >= freq_vector[i] & f <= freq_vector[i+1]
    H[idx] <- amp_vector[i]
  }
  
  # Get impulse response
  h <- fft(H, inverse = TRUE) / n_fft
  h <- Re(h)
  
  # Apply Hamming window
  w <- hamming(n_fft)
  h <- h * w
  
  return(h)
}
# Function to plot frequency response with -3dB line
plot_filter_freq_response <- function(h, fs) {
  # Calculate frequency response using FFT
  n_fft <- 64
  H <- fft(c(h, rep(0, n_fft - length(h))))
  freq <- seq(0, fs/2, length.out = floor(n_fft/2) + 1)
  mag <- abs(H[1:length(freq)])
  
  # Create plot data
  df <- tibble(
    frequency = freq,
    magnitude = mag/max(mag)  # Normalize magnitude
  ) 
  
  # Create plot
  ggplot(df, aes(x = frequency, y = magnitude)) +
    # Add -3dB line
    geom_hline(yintercept = 1/sqrt(2), linetype = "dashed", color = "red", alpha = 0.5) +
    geom_line() +
    # Add text label for -3dB line
    annotate("text", x = max(freq) * 0.8, y = 1/sqrt(2) + 0.05, 
             label = "-3dB", color = "red") +
    labs(x = "Frequency (Hz)", 
         y = "Magnitude") +
    ylim(0, 1.1) +
    theme_minimal()
}

process_rt_signal_eeglab <- function(rt_data, fs = 100) {
  # Function to process RT data following EEGLAB's methodology
  
  # 1. Zero padding (50 points pre and post)
  pad_length <- 50
  padded_data <- c(rep(0, pad_length), rt_data, rep(0, pad_length))
  
  # 2. Apply Hanning window
  window <- hanning(length(padded_data))
  windowed_data <- padded_data * window
  
  # 3. Create EEGLAB-style filters for each band
  filters <- list(
    band1 = design_eeglab_fir(0, 2, fs, order = 10),    # 0-2 Hz
    band2 = design_eeglab_fir(2, 5, fs, order = 10),    # 2-5 Hz
    band3 = design_eeglab_fir(8, 20, fs, order = 10)    # 8-20 Hz
  )
  
  # 4. Two-pass filtering for each band
  filter_two_pass <- function(data, filter_coef) {
    # Forward pass
    forward <- stats::filter(data, filter_coef)
    # Reverse the data
    reverse <- rev(forward)
    # Backward pass
    backward <- stats::filter(reverse, filter_coef)
    # Reverse again to get original orientation
    return(rev(backward))
  }
  
  # Apply filters to windowed data
  filtered_bands <- lapply(filters, function(filter) {
    filter_two_pass(windowed_data, filter)
  })
  
  # Remove padding from filtered data
  filtered_bands <- lapply(filtered_bands, function(band) {
    band[(pad_length + 1):(length(band) - pad_length)]
  })
  
  # Return results
  list(
    original = rt_data,
    filtered_bands = filtered_bands,
    filters = filters
  )
}

# Example usage
set.seed(123)
fs <- 60  # 60 Hz sampling rate
t <- seq(0, 1, 1/fs)  # 1 second of data

# Create synthetic RT data
rt_data <- 0.5 * sin(2*pi*1*t) +    # 1 Hz component
  0.3 * sin(2*pi*3.5*t) +    # 3.5 Hz component
  0.2 * sin(2*pi*18*t)       # 18 Hz component

# Process the data
results <- process_rt_signal_eeglab(rt_data, fs)  # Note: removed padding here

# Create plots
p1 <- plot_filter_freq_response(results$filters$band1, fs) + 
  ggtitle("0-2 Hz Band Filter")
p2 <- plot_filter_freq_response(results$filters$band2, fs) + 
  ggtitle("2-5 Hz Band Filter")
p3 <- plot_filter_freq_response(results$filters$band3, fs) + 
  ggtitle("8-20 Hz Band Filter")

# Plot filtered signals - ensure all vectors have the same length
signal_df <- tibble(
  time = t,
  original = results$original,
  band1 = results$filtered_bands$band1,
  band2 = results$filtered_bands$band2,
  band3 = results$filtered_bands$band3
)

# Create separate plots for each band's time domain response
p4 <- ggplot(signal_df, aes(x = time)) +
  geom_line(aes(y = original), alpha = 0.3) +
  geom_line(aes(y = band1), color = "blue") +
  labs(title = "0-2 Hz Band",
       x = "Time (s)",
       y = "Amplitude") +
  theme_minimal()

p5 <- ggplot(signal_df, aes(x = time)) +
  geom_line(aes(y = original), alpha = 0.3) +
  geom_line(aes(y = band2), color = "red") +
  labs(title = "2-5 Hz Band",
       x = "Time (s)",
       y = "Amplitude") +
  theme_minimal()

p6 <- ggplot(signal_df, aes(x = time)) +
  geom_line(aes(y = original), alpha = 0.3) +
  geom_line(aes(y = band3), color = "green") +
  labs(title = "8-20 Hz Band",
       x = "Time (s)",
       y = "Amplitude") +
  theme_minimal()

# Arrange plots in a 3x2 grid
(p1 + p4) / (p2 + p5) / (p3 + p6)

# Print filter coefficients
cat("\nFilter coefficients:\n")
cat("\n0-2 Hz band:\n")
print(round(results$filters$band1, 6))
cat("\n2-5 Hz band:\n")
print(round(results$filters$band2, 6))
cat("\n8-20 Hz band:\n")
print(round(results$filters$band3, 6))