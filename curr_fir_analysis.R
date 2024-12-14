library(tidyverse)
library(signal)
library(pracma)
library(patchwork)

# Post-2015 EEGLAB FIR design (firls approach)
design_eeglab_fir_post2015 <- function(low_freq, high_freq, fs, force_filtorder = NULL) {
  nyquist <- fs/2
  if (is.null(low_freq)) low_freq <- 0
  if (is.null(high_freq)) high_freq <- nyquist
  
  if (!is.null(force_filtorder)) {
    filtorder <- force_filtorder
  } else {
    if (low_freq > 0) {
      filtorder <- ceiling(3.3 * fs / low_freq)
    } else {
      filtorder <- ceiling(3.3 * fs / high_freq)
    }
  }
  
  if (filtorder %% 2 == 0) filtorder <- filtorder + 1
  
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
  
  h <- signal::fir2(filtorder, freq_vector/nyquist, amp_vector)
  
  return(list(
    h = h,
    freq_vector = freq_vector,
    amp_vector = amp_vector,
    filtorder = filtorder
  ))
}

plot_filter_freq_response <- function(filter_result, fs, 
                                      show_x_lab = FALSE,
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
  
  db_threshold <- 1/sqrt(2)
  crossings <- which(diff(df$magnitude > db_threshold) != 0)
  
  db_freqs <- numeric()
  for(i in crossings) {
    x1 <- df$frequency[i]
    x2 <- df$frequency[i + 1]
    y1 <- df$magnitude[i]
    y2 <- df$magnitude[i + 1]
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
  
  # Create custom breaks for x-axis that include both standard breaks and -3dB points
  standard_breaks <- seq(0, 25, 5)
  all_breaks <- sort(unique(c(standard_breaks, db_freqs)))
  
  # Create custom labels with colors
  break_labels <- all_breaks
  
  p <- ggplot() +
    theme_bw() +
    theme(
      text = element_text(family = "Times New Roman", size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(color = ifelse(all_breaks %in% db_freqs, "darkred", "black")),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      panel.border = element_rect(color = "black", linewidth = 1),
      plot.title = element_blank(),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    ) +
    geom_hline(yintercept = db_threshold, linetype = "dashed", 
               color = "darkgray", linewidth = 0.5) +
    geom_line(data = df, aes(x = frequency, y = magnitude), linewidth = 0.7) +
    geom_point(data = design_points, aes(x = frequency, y = magnitude), 
               color = "black", size = 2, shape = 21, fill = "white") +
    geom_point(data = db_points, aes(x = frequency, y = magnitude),
               color = "darkred", size = 2) +
    geom_segment(data = db_points,
                 aes(x = frequency, xend = frequency,
                     y = -0.02, yend = magnitude),
                 color = "darkred", linetype = "dotted", linewidth = 0.5) +
    annotate("text", x = 1, y = db_threshold, 
             label = "-3dB", hjust = -0.2, vjust = -0.5,
             size = 3, color = "darkgray") +
    labs(x = if(show_x_lab) "Frequency (Hz)" else "",
         y = paste0("Magnitude\n", y_text)) +
    scale_y_continuous(limits = c(-0.1, 1.1), 
                       breaks = seq(0, 1, 0.2),
                       expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 25),
                       breaks = all_breaks,
                       labels = sprintf("%.1f", break_labels),
                       expand = c(0.02, 0))
  
  return(p)
}

# Set parameters
fs <- 60

# Design filters (only least-squares)
filters_post2015 <- list(
  band1 = design_eeglab_fir_post2015(0, 2, fs, force_filtorder = 11),
  band2 = design_eeglab_fir_post2015(2, 5, fs, force_filtorder = 11),
  band3 = design_eeglab_fir_post2015(8, 20, fs, force_filtorder = 11)
)

# Create plots
p1 <- plot_filter_freq_response(filters_post2015$band1, fs, FALSE, "0-2 Hz")
p2 <- plot_filter_freq_response(filters_post2015$band2, fs, FALSE, "2-5 Hz")
p3 <- plot_filter_freq_response(filters_post2015$band3, fs, TRUE, "8-20 Hz")

# Create vertical layout
combined_plot <- p1 / p2 / p3 +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Least-Squares FIR Filter Design",
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