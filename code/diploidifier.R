# Haploid combiner
#
# Input: VCF with two subjects' X chromosome genotypes
# Output: VCF with pseudo-diploid constructed using haploids from input
#
# Example:
# Rscript diploidifier.R input.vcf.gz output.vcf.gz

library(tidyverse)
library(vcfR)

args <- commandArgs(trailingOnly = T)
input_file <- args[1]
output_file <- args[2]

#' Create shuffled diploid from two genotypes
#'
#' @param g1 A genotype (0 or 1)
#' @param g2 A genotype (0 or 1)
#' @return An unphased diploid in VCF format (1/0, 1/1, 0/1, or 0/0)
unphaser <- function(g1, g2){
  coin_flip <- ifelse(runif(1)>0.5, TRUE, FALSE)
  if(coin_flip){
    val <- paste0(g1, "/", g2)
  } else{
    val <- paste0(g2, "/", g1)
  }
  return(val)
}

vcf <- read.vcfR(input_file, verbose = FALSE)
vcf_truth <- vcf

df <- vcf@gt %>% as.data.frame()
s1 <- names(df)[2]
s2 <- names(df)[3]
df <- df %>%
  rowwise() %>%
  mutate(FAKE001 = unphaser(str_sub( !!as.name(s1), 1, 1), str_sub( !!as.name(s2), 1, 1)) ) %>%
  ungroup() %>%
  select(FORMAT, FAKE001)

vcf@gt <- as.matrix(df)
write.vcf(vcf, paste0(output_file, "_test.vcf.gz"))

df <- vcf_truth@gt %>% as.data.frame()
s1 <- names(df)[2]
s2 <- names(df)[3]
df <- df %>%
  rowwise() %>%
  mutate(FAKE001 = paste0(str_sub( !!as.name(s1), 1, 1), "|", str_sub( !!as.name(s2), 1, 1))) %>%
  ungroup() %>%
  select(FORMAT, FAKE001)

vcf_truth@gt <- as.matrix(df)
write.vcf(vcf_truth, paste0(output_file, "_truth.vcf.gz"))
