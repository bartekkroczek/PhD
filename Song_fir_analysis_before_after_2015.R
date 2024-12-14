library(tidyverse)
library(signal)
library(pracma)
library(patchwork)

# Pre-2015 EEGLAB FIR design (windowed sinc)
design_eeglab_fir_pre2015 <- function(low_freq, high_freq, fs, revfilt = 0, force_filtorder = NULL) {
  nyquist <- fs/2
  if (is.null(low_freq)) low_freq <- 0
  if (is.null(high_freq)) high_freq <- nyquist
  
  # Use forced filter order if provided, otherwise calculate as before
  if (!is.null(force_filtorder)) {
    filtorder <- force_filtorder
  } else {
    if (low_freq > 0) {
      filtorder <- 3 * floor(fs / low_freq)
    } else {
      filtorder <- 3 * floor(fs / high_freq)
    }
  }
  
  if (filtorder %% 2 == 0) filtorder <- filtorder + 1
  if (low_freq > 0) {
    trans_band_low <- if (low_freq <= 2) {
      min(2, low_freq * 0.5)
    } else if (low_freq <= 5) {
      min(3, low_freq * 0.35)
    } else if (low_freq <= 10) {
      min(4, low_freq * 0.25)
    } else {
      low_freq * 0.2
    }
  }
  
  if (high_freq < nyquist) {
    trans_band_high <- if (high_freq <= 2) {
      min(2, high_freq * 0.5)
    } else if (high_freq <= 5) {
      min(3, high_freq * 0.35)
    } else if (high_freq <= 10) {
      min(4, high_freq * 0.25)
    } else {
      high_freq * 0.2
    }
  }
  
  if (low_freq <= 0) {
    freq_vector <- c(0, 
                     high_freq - trans_band_high,
                     high_freq,
                     high_freq + trans_band_high,
                     nyquist)
    amp_vector <- c(1, 1, 0.5, 0, 0)
  } else if (high_freq >= nyquist) {
    freq_vector <- c(0,
                     low_freq - trans_band_low,
                     low_freq,
                     low_freq + trans_band_low,
                     nyquist)
    amp_vector <- c(0, 0, 0.5, 1, 1)
  } else {
    freq_vector <- c(0,
                     low_freq - trans_band_low,
                     low_freq,
                     low_freq + trans_band_low,
                     high_freq - trans_band_high,
                     high_freq,
                     high_freq + trans_band_high,
                     nyquist)
    amp_vector <- c(0, 0, 0.5, 1, 1, 0.5, 0, 0)
  }
  
  freq_vector <- freq_vector / nyquist
  n_fft <- filtorder + 1
  f <- seq(0, 1, length.out = n_fft)
  H <- pracma::interp1(freq_vector, amp_vector, f, method = "linear")
  H[is.na(H)] <- 0
  
  h <- fft(H, inverse = TRUE) / n_fft
  h <- Re(h)
  
  n <- length(h)
  w <- 0.54 - 0.46 * cos(2 * pi * (0:(n-1)) / (n-1))
  h <- h * w
  
  if (revfilt == 1) amp_vector <- 1 - amp_vector
  
  return(list(
    h = h,
    freq_vector = freq_vector * nyquist,
    amp_vector = amp_vector,
    filtorder = filtorder
  ))
}

# Post-2015 EEGLAB FIR design (firls approach)
design_eeglab_fir_post2015 <- function(low_freq, high_freq, fs, force_filtorder = NULL) {
  nyquist <- fs/2
  if (is.null(low_freq)) low_freq <- 0
  if (is.null(high_freq)) high_freq <- nyquist
  
  # Use forced filter order if provided, otherwise calculate as before
  if (!is.null(force_filtorder)) {
    filtorder <- force_filtorder
  } else {
    if (low_freq > 0) {
      filtorder <- ceiling(3.3 * fs / low_freq)
    } else {
      filtorder <- ceiling(3.3 * fs / high_freq)
    }
  }
  
  # Ensure odd filter order
  if (filtorder %% 2 == 0) filtorder <- filtorder + 1
  
  # Transition bands (firfilt uses different transition widths)
  if (low_freq <= 0) {  # Lowpass
    freq_vector <- c(0, 
                     0.85 * high_freq,
                     high_freq,
                     nyquist)
    amp_vector <- c(1, 1, 0, 0)
  } else if (high_freq >= nyquist) {  # Highpass
    freq_vector <- c(0,
                     low_freq,
                     1.15 * low_freq,
                     nyquist)
    amp_vector <- c(0, 0, 1, 1)
  } else {  # Bandpass
    freq_vector <- c(0,
                     low_freq,
                     1.15 * low_freq,
                     0.85 * high_freq,
                     high_freq,
                     nyquist)
    amp_vector <- c(0, 0, 1, 1, 0, 0)
  }
  
  # Design filter using least squares
  h <- signal::fir2(filtorder, freq_vector/nyquist, amp_vector)
  
  return(list(
    h = h,
    freq_vector = freq_vector,
    amp_vector = amp_vector,
    filtorder = filtorder
  ))
}

# Enhanced plotting function with -3dB markers on x-axis
plot_filter_freq_response <- function(filter_result, fs, 
                                      show_x_lab = FALSE, 
                                      show_y_lab = FALSE,
                                      y_text = "") {
  n_fft <- 512
  h <- filter_result$h
  H <- fft(c(h, rep(0, n_fft - length(h))))
  freq <- seq(0, fs/2, length.out = floor(n_fft/2) + 1)
  mag <- abs(H[1:length(freq)])
  
  df <- tibble(
    frequency = freq,
    magnitude = mag/max(mag)
  )
  
  # Find -3dB points (magnitude ≈ 0.707)
  db_threshold <- 1/sqrt(2)
  crossings <- which(diff(df$magnitude > db_threshold) != 0)
  
  # Linear interpolation to find exact frequencies
  db_freqs <- numeric()
  for(i in crossings) {
    x1 <- df$frequency[i]
    x2 <- df$frequency[i + 1]
    y1 <- df$magnitude[i]
    y2 <- df$magnitude[i + 1]
    
    # Linear interpolation
    x_intercept <- x1 + (db_threshold - y1) * (x2 - x1)/(y2 - y1)
    db_freqs <- c(db_freqs, x_intercept)
  }
  
  design_points <- tibble(
    frequency = filter_result$freq_vector,
    magnitude = filter_result$amp_vector
  )
  
  db_points <- tibble(
    frequency = db_freqs,
    magnitude = rep(db_threshold, length(db_freqs))
  )
  
  p <- ggplot() +
    theme_bw() +
    theme(
      text = element_text(family = "Times New Roman", size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      panel.border = element_rect(color = "black", linewidth = 1),
      plot.title = element_blank(),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2),
      axis.title.x = if(show_x_lab) element_text() else element_blank(),
      axis.title.y = if(show_y_lab) element_text() else element_blank()
    ) +
    # Add -3dB line
    geom_hline(yintercept = db_threshold, linetype = "dashed", 
               color = "darkgray", linewidth = 0.5) +
    geom_line(data = df, aes(x = frequency, y = magnitude), linewidth = 0.7) +
    geom_point(data = design_points, aes(x = frequency, y = magnitude), 
               color = "black", size = 2, shape = 21, fill = "white") +
    # Add -3dB points
    geom_point(data = db_points, aes(x = frequency, y = magnitude),
               color = "darkred", size = 2) +
    # Add vertical lines to x-axis for -3dB points
    geom_segment(data = db_points,
                 aes(x = frequency, xend = frequency,
                     y = -0.02, yend = magnitude),
                 color = "darkred", linetype = "dotted", linewidth = 0.5) +
    # Add -3dB label on y-axis
    annotate("text", x = 1, y = db_threshold, 
             label = "-3dB", hjust = -0.2, vjust = -0.5,
             size = 3, color = "darkgray") +
    # Add frequency labels on x-axis
    geom_text(data = db_points,
              aes(x = frequency, y = -0.05,
                  label = sprintf("%.1f", frequency)),
              size = 2.5, color = "darkred", vjust = 1) +
    labs(x = if(show_x_lab) "Frequency (Hz)" else "",
         y = if(show_y_lab) paste0("Magnitude\n", y_text) else "") +
    scale_y_continuous(limits = c(-0.1, 1.1), 
                       breaks = seq(0, 1, 0.2),
                       expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 25),
                       breaks = seq(0, 25, 5),
                       expand = c(0.02, 0))
  
  return(p)
}

# Set parameters
fs <- 60

# Design filters
filters_pre2015 <- list(
  band1 = design_eeglab_fir_pre2015(0, 2, fs, force_filtorder = 10),
  band2 = design_eeglab_fir_pre2015(2, 5, fs, force_filtorder = 10),
  band3 = design_eeglab_fir_pre2015(8, 20, fs, force_filtorder = 10)
)

filters_post2015 <- list(
  band1 = design_eeglab_fir_post2015(0, 2, fs, force_filtorder = 10),
  band2 = design_eeglab_fir_post2015(2, 5, fs, force_filtorder = 10),
  band3 = design_eeglab_fir_post2015(8, 20, fs, force_filtorder = 10)
)

# Create plots with band labels in y-axis
p1_pre <- plot_filter_freq_response(filters_pre2015$band1, fs, FALSE, TRUE, "0-2 Hz")
p2_pre <- plot_filter_freq_response(filters_pre2015$band2, fs, FALSE, TRUE, "2-5 Hz")
p3_pre <- plot_filter_freq_response(filters_pre2015$band3, fs, TRUE, TRUE, "8-20 Hz")

p1_post <- plot_filter_freq_response(filters_post2015$band1, fs, FALSE, FALSE)
p2_post <- plot_filter_freq_response(filters_post2015$band2, fs, FALSE, FALSE)
p3_post <- plot_filter_freq_response(filters_post2015$band3, fs, TRUE, FALSE)

# Create compact layout
combined_plot <- (p1_pre | p1_post) /
  (p2_pre | p2_post) /
  (p3_pre | p3_post) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Windowed-Sinc Design                                    Least-Squares Design",
    theme = theme(
      plot.title = element_text(
        family = "Times New Roman",
        size = 12,
        hjust = 0.5,
        margin = margin(t = 5, b = 5)
      ),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  )

# Display the plot
combined_plot