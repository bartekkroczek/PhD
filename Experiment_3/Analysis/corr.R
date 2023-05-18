library(tidyverse)
library(readr)

list_of_files <- list.files(path = '../Data/saccade_task/csv',
                            pattern = '*.csv$',
                            full.names = TRUE)

df <- readr::read_csv(list_of_files)

