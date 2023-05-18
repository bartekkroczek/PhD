library(tidyverse)
library(readr)

list_of_files <- list.files(path = '../Data/saccade_task/',
                            pattern = '*.csv$',
                            full.names = TRUE)

df <- readr::read_csv(list_of_files,
                      col_types = cols(.default = col_character(),
                                       Block_no = col_integer(),
                                       Trial_no = col_integer(),
                                       Block_type = col_factor(),
                                       Trial_type = col_factor(),
                                       CSI = col_integer(),
                                       Rt = col_double(),
                                       Corr = col_logical(),
                                       "Stimulus Time" = col_integer()))

