# Annotate CpG switches with methylation data

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
id <- args[2]

library(tidyverse)

f_name <- paste0(input_dir, "/switch_", id, ".csv" )
out_name <- paste0(input_dir, "/switch_meth_", id, ".csv")
df <- read_csv(f_name)

df <- df %>%
  mutate(center_nuc = str_sub(motif_start, 2, 2)) %>%
  mutate(start = ifelse(center_nuc == "C", pos_start - 1, ifelse(center_nuc == "G", pos_start - 2, NA)))

meth_df <- read_csv("/net/snowwhite/home/beckandy/research/phasing/output/meth_locations.csv")

df <- left_join(df, {meth_df %>% select(start, meth)}, by = "start")

write_csv(df, out_name)
