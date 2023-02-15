# 16 November 2022
# The purpose of this code is to annotate the list of heterozygous sites in each pseudodiploid
# with whether or not the position is a CpG

library(tidyverse)

input_dir <- "/net/snowwhite/home/beckandy/research/phasing/output/2023_switch_errors/het_loc/"
out_dir <- paste0(input_dir, "annotated/")

df_cpg <- read_tsv("/net/snowwhite/home/beckandy/research/phasing/data/cpg_pos_mask.bed",
                   col_names = c("chrom", "start", "end"))

load_pair <- function(pair_id, input_dir ="/net/snowwhite/home/beckandy/research/phasing/output/2023_switch_errors/het_loc/"){
  f_name <- paste0(input_dir, "pair_", pair_id, "_het_loc.txt")
  df <- read_tsv(f_name, col_names = c("chrom", "pos", "gt")) %>%
    select(chrom, pos, gt)
  return(df)
}
## work out a single example
df <- load_pair(1)
df$cpg <- as.numeric(df$pos %in% df_cpg$start)

## wow, that was hard :p
## do it for all and save

append_cpg <- function(pair_id, input_dir ="/net/snowwhite/home/beckandy/research/phasing/output/2023_switch_errors/het_loc/"){
  df <- load_pair(pair_id, input_dir)
  df$cpg <- as.numeric(df$pos %in% df_cpg$start)
  return(df)
}

for(i in 1:400){
  df <- append_cpg(i)
  f_name <- paste0(out_dir, "pair_", i, "_het_loc.txt")
  write_tsv(df, f_name)
}
