library(ggplot2)
library(jtools)

# Set up PDI range (50-3000 ms)
pdi_range <- seq(50, 3000, by = 25)

# Define theoretical models
# 1. Decay hypothesis: linear decline
decay_model <- function(pdi, start_acc = 85, slope = -0.003) {
  start_acc + slope * pdi
}

# 2. Preparatory enhancement: linear increase
prep_model <- function(pdi, start_acc = 85, slope = 0.002) {
  start_acc + slope * pdi
}

# 3. Oscillatory attention: damped sine wave
osc_model <- function(pdi, A_avg = 85, A_amp = 8, freq = 0.004, decay_rate = 0.0003) {
  A_avg + A_amp * sin(2 * pi * freq * pdi) * exp(-decay_rate * pdi)
}

# 4. Null hypothesis: flat line
null_model <- function(pdi, mean_acc = 85) {
  rep(mean_acc, length(pdi))
}

# Create theoretical curves dataframe
theories_df <- data.frame(
  PDI = rep(pdi_range, 4),
  Accuracy = c(
    decay_model(pdi_range),
    prep_model(pdi_range),
    osc_model(pdi_range),
    null_model(pdi_range)
  ),
  Hypothesis = rep(c("Deterioration over time", "Improvement over time", "Behavioral Oscillation", "No Effect"), 
                   each = length(pdi_range))
)

# Define colors
colors <- c("Deterioration over time" = "#D55E00", 
            "Improvement over time" = "#009E73", 
            "Behavioral Oscillation" = "#0072B2", 
            "No Effect" = "#999999")

# Create plot
ggplot(theories_df, aes(x = PDI, y = Accuracy, color = Hypothesis)) +
  geom_line(size = 1.2, alpha = 0.9) +
  
  scale_color_manual(values = colors) +
  
  labs(
    x = "Preparatory Delay Interval (ms)",
    y = "Mean Discrimination Accuracy (%)",
    color = "Hypothesis"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  ) +
  coord_cartesian(ylim = c(75, 95))








################ second plot 


library(ggplot2)
library(jtools)

# Two curves starting at the same point (x = 0)
set.seed(1)
x <- seq(0, 3000, by = 5)
base <- 0.55  # common starting accuracy

# Line 1: linear increase over time
y_inc <- base + (0.92 - base) * (x / 3000)

# Line 2: fast rise (pre-600 ms), then flat (plateau) after 600 ms
tau <- 150  # controls the speed of early rise
y600 <- base + 0.35 * (1 - exp(-600 / tau))
y_ab <- ifelse(x <= 600, base + 0.35 * (1 - exp(-x / tau)), y600)

df_inc <- data.frame(x = x, y = y_inc, model = "Increment over time hypothesis")
df_ab  <- data.frame(x = x, y = y_ab,  model = "Attentional blink hypothesis")
df <- rbind(df_inc, df_ab)

ggplot(df, aes(x = x, y = y, color = model, linetype = model)) +
  # Dimmed zone for 0–600 ms
  geom_rect(aes(xmin = 0, xmax = 600, ymin = -Inf, ymax = Inf),
            fill = "grey30", alpha = 1, inherit.aes = FALSE) +
  annotate("text", x = 300, y = 0.96, label = "Attentional blink\n effects zone",
           color = "black", size = 5) +
  # Add a tiny right-side margin while keeping axis labeled 0–3000
  scale_x_continuous(limits = c(0, 3020), breaks = c(0, 1000, 2000, 3000), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 1.0), breaks = NULL, labels = NULL, expand = c(0, 0)) +
  scale_color_manual(values = c("Increment over time hypothesis" = "#009E73",
                                "Attentional blink hypothesis" = "#0072B2")) +
  scale_linetype_manual(values = c("Increment over time hypothesis" = "solid",
                                   "Attentional blink hypothesis" = "solid")) +
  labs(x = "Preparatory delay interval (ms)",
       y = "Attention accuracy",
       color = NULL,
       linetype = NULL) +
  theme_apa() +
  theme(
    legend.position = "top",
    panel.border = element_blank(),                  # remove full panel border
    axis.line = element_line(color = "black"),       # draw axes
    axis.line.x.top = element_blank(),               # remove top axis line
    axis.line.y.right = element_blank(),             # remove right axis line
    plot.margin = margin(8, 18, 8, 8)                # add a bit of outer margin (t, r, b, l)
  ) +
  coord_cartesian(xlim = c(0, 3000), clip = "off")   # keep labels 0–3000 and avoid clipping


#####



