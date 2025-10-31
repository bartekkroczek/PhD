library(tidyverse)

vec <- read_csv("cleaned_cognitive_data.csv")$PART_ID

# Extract letter (K or M) - 3rd character from end
sex <- str_sub(vec, -3, -3)

# Extract last two digits as numbers
age <- as.numeric(str_sub(vec, -2, -1))

# Results
print(table(sex))
print(table(sex)/length(vec))
print(summary(age))
print(sd(age))