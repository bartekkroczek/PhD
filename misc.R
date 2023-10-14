library(tidyverse)
df <- read_csv("Experiment_1/Data/127K22_795_beh.csv")
view(df)

res <- df |> 
  filter(Trial_type == 'experiment') 

unique(res$CSI)