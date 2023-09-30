# Description of the script

library(optparse)
library(tidyverse)
library(openxlsx)

combined.dd <- read.xlsx("/home/fdelogu/biowide_net/results/split_nets/combined_nets.xlsx", colNames = T)

combined.dd %>%
  filter(tech!="sparcc") %>%
  ggplot(aes(x=p.value, y=stat, color=tech)) +
  geom_point(alpha=0.1)