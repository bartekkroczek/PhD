library(tidyverse)
library(corrr)
library(readr)
library(stats)
library(psych)
library(ggstatsplot)
library(cowplot)
library(ggpubr)
library(jtools)



list_of_files <- list.files(path = '../Data/saccade_task/csv',
                            pattern = ".*csv$",
                            full.names = TRUE)

df <- readr::read_csv(list_of_files)


ms_per_frame <- 1000. / 60

adapt_per_task <- df |>
  filter(str_detect(PART_ID, "PB"), Trial_type == "exp") |>
  group_by(PART_ID, Block_no, Block_type) |>
  summarise(adapt_val = mean(tail(`Stimulus Time`, 6)) * ms_per_frame) |>
  group_by(PART_ID, Block_type) |>
  summarise(blk_adapt_val = mean(adapt_val)) |>
  pivot_wider(names_from = Block_type,
              values_from = blk_adapt_val) |>
  mutate(DIFF = AS - PS) |>
  filter(AS < 500, PS < 250, DIFF > 10, DIFF < 250) 
#   ungroup() |>
# select(AS, PS, DIFF) |>
# pairs.panels(gap = 0, pch = 21, lm = TRUE, ci = TRUE)
# ggplot(mapping = aes(x=AS, y=DIFF)) +
# geom_point() +
# geom_smooth(method = "lm")


# Kogo wyrzucam
#   PART_ID     AS    PS
# 1 086M34    301.  288.
# 2 PB014M27  637.  135.
# 3 PB026K36  331.  269.
# 4 PB060K22  612.  349.

as_ps <- ggplot(adapt_per_task, mapping = aes(x = PS, y = AS)) +
  geom_point(size = 2.5, color = "blue") +
  geom_smooth(method = "lm", color = "black") +
  xlim(c(100, 240)) +
  ylim(c(150, 410)) +
  labs(x = "Prosaccade duration", y = "Antisaccade duration") +
  stat_cor(
    method = "pearson",
    label.x = 110,
    label.y = 350,
    p.accuracy = 0.001,
    size = 6
  ) +
  theme_apa(x.font.size = 18, y.font.size = 18)

ps_diff <- ggplot(adapt_per_task, mapping = aes(x = PS, y = DIFF)) +
  geom_point(size = 2.5, color = "blue") +
  geom_smooth(method = "lm", color = "black") +
  xlim(c(100, 240)) +
  ylim(c(10, 200)) +
  labs(x = "Prosaccade duration", y = "Antisaccade − Prosaccade") +
  stat_cor(
    method = "pearson",
    label.x = 105,
    label.y = 160,
    p.accuracy = 0.001,
    size = 6
  ) +
  theme_apa(x.font.size = 18, y.font.size = 18)

as_diff <- ggplot(adapt_per_task, mapping = aes(x = AS, y = DIFF)) +
  geom_point(size = 2.5, color = "blue") +
  geom_smooth(method = "lm", color = "black") +
  xlim(c(150, 410)) +
  ylim(c(10, 200)) +
  labs(x = "Antisaccade duration", y = "Antisaccade − Prosaccade") +
  stat_cor(
    method = "pearson",
    label.x = 160,
    label.y = 160,
    p.accuracy = 0.001,
    size = 6
  ) +
  theme_apa(x.font.size = 18, y.font.size = 18)

plot_grid(as_ps, ps_diff, as_diff, nrow = 1, labels = c('A', 'B', 'C'))

ggscatterstats(
  data = adapt_per_task, 
  x = PS, 
  y = AS,
  xlab = "Prosaccade Duration",
  ylab = "Antisaccade Duration")



par(
  pty = 's',
  mfrow = c(1, 3),
  mar = c(4, 4, 1, 1),
  mgp = c(2, 1, 0),
  cex.lab = 1.2
)

plot(
  allTimes[, 1:2],
  typ = 'n',
  ylim = c(0, 200),
  xlim = c(0, 200),
  xlab = "Prosaccade Duration",
  ylab = "Antisaccade Duration"
)
segments(m[, 1, 1], allTimes[, 2], m[, 2, 1], allTimes[, 2], col = 'darkblue')
segments(allTimes[, 1], m[, 1, 2], allTimes[, 1], m[, 2, 2], col = 'darkblue')
points(allTimes[, 1:2],
       pch = 19,
       col = rgb(.5, .5, 1, .3),
       cex = 2)
abline(lm(allTimes[, 2] ~ allTimes[, 1]), lty = 2)
abline(0, 1)
text(100, 50, bquote(R == .(cors[1, 2])), cex = 1.2)
mtext(side = 3,
      adj = .5,
      cex = 1.2,
      "A.")

plot(
  allTimes[, c(1, 3)],
  pch = 19,
  col = 'blue',
  ylim = c(0, 120),
  xlim = c(0, 120),
  xlab = "Prosaccade Duration",
  ylab = "Antisaccade - Prosaccade"
)
abline(lm(allTimes[, 3] ~ allTimes[, 1]), lty = 2)
abline(0, 1)
text(100, 50, bquote(R == .(cors[1, 3])), cex = 1.2)
mtext(side = 3,
      adj = .5,
      cex = 1.2,
      "B.")

plot(
  allTimes[, c(2, 3)],
  pch = 19,
  col = 'blue',
  ylim = c(0, 200),
  xlim = c(0, 200),
  xlab = "Antisaccade Duration",
  ylab = "Antisaccade - Prosaccade"
)
abline(lm(allTimes[, 3] ~ allTimes[, 2]), lty = 2)
abline(0, 1)
text(50, 100, bquote(R == .(cors[2, 3])), cex = 1.2)
mtext(side = 3,
      adj = .5,
      cex = 1.2,
      "C.")