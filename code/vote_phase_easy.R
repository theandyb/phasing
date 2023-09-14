library(tidyverse)
library(vcfR)

args <- commandArgs(trailingOnly = T)
id <- args[1]
output_file <- args[2]

vcf_dir <- "/net/snowwhite/home/beckandy/research/phasing/output/final_switch_errors/vote/easy/"
subj_name <- paste0("FAKE", id)

# read in phased files
vcf_eagle <- read.vcfR(paste0(vcf_dir, "eagle_", id, ".vcf.gz"),
                       verbose = FALSE )
vcf_beagle <- read.vcfR(paste0(vcf_dir, "beagle_", id, ".vcf.gz"),
                        verbose = FALSE )
vcf_shapeit <- read.vcfR(paste0(vcf_dir, "shapeit_", id, ".vcf.gz"),
                         verbose = FALSE )

# get haploids as data.frame objects
df_eagle <- vcf_eagle@gt %>% as.data.frame()
df_beagle <- vcf_beagle@gt %>% as.data.frame()
df_shapeit <- vcf_shapeit@gt %>% as.data.frame()

#let's split the GT calls into two columns
df_eagle <- df_eagle %>%
  separate(!!subj_name, c("v1", "v2"), sep = "\\|")
df_beagle <- df_beagle %>%
  separate(!!subj_name, c("v1", "v2"), sep = "\\|")
df_shapeit <- df_shapeit %>%
  separate(!!subj_name, c("v1", "v2"), sep = "\\|")

# template for new haploids
v1 <- rep(0, length(df_shapeit$v1))
v2 <- rep(0, length(df_shapeit$v2))

for(i in 1:length(df_shapeit$v1)){
  # first check for agreement
  if(df_eagle$v1[i] == df_beagle$v1[i] & df_eagle$v1[i] == df_shapeit$v1[i]){
    # match -> no vote necessary
    v1[i] <- df_shapeit$v1[i]
    v2[i] <- df_shapeit$v2[i]
  } else {
    # find vote winner at position
    if(df_eagle$v1[i] == df_shapeit$v1[i]){
      v1[i] <- df_shapeit$v1[i]
      v2[i] <- df_shapeit$v2[i]
    } else if (df_beagle$v1[i] == df_shapeit$v1[i]){
      v1[i] <- df_shapeit$v1[i]
      v2[i] <- df_shapeit$v2[i]
    } else{
      v1[i] <- df_eagle$v1[i]
      v2[i] <- df_eagle$v2[i]
    }
  }
}

vcf_new <- vcf_eagle
df_new <- vcf_eagle@gt %>% as.data.frame()
df_new[ , subj_name] <- paste0(v1, "|", v2)
vcf_new@gt <- as.matrix(df_new)
write.vcf(vcf_new, output_file)
