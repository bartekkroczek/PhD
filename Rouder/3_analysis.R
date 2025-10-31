
  
  ## Results 
  
  
cors <- tmp_df |>
  select(pro, anti, diff) |>
  cor() |>
  round(2)

create_saccade_plots(tmp_df, cors, contrast = "diff", y_lab = "Antisaccade - Prosaccade")


# Calculate correlations

cors <- tmp_df |>
  select(pro, anti, processing_speed) |>
  cor() |>
  round(2)

create_saccade_plots(tmp_df, cors, contrast = "processing_speed", y_lab = "Processing Speed")





# Is prosaccade in fact a processing speed?


ggplot(tmp_df, aes(x = processing_speed, y = pro)) +
  geom_point(color = "blue", size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linetype = 2, color = "black") +
  annotate("text", x = mean(tmp_df$processing_speed), y = max(tmp_df$pro) * 0.9, label = paste0("R = ", round(cor(tmp_df$processing_speed, tmp_df$pro), 2))) +
  labs(x = "Processing Speed", y = "Prosacade") +
  #  coord_cartesian(xlim = c(100, 450), ylim = c(0, 250)) +
  theme_minimal()

# Is cognitive control in fact a processing speed?


ggplot(tmp_df, aes(x = processing_speed, y = diff)) +
  geom_point(color = "blue", size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linetype = 2, color = "black") +
  annotate("text", x = mean(tmp_df$processing_speed), y = max(tmp_df$diff) * 0.9, label = paste0("R = ", round(cor(tmp_df$processing_speed, tmp_df$diff), 2))) +
  labs(x = "Processing Speed", y = "Antisaccade - Prosacade") +
  #  coord_cartesian(xlim = c(100, 450), ylim = c(0, 250)) +
  theme_minimal()



