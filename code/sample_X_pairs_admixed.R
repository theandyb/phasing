# Sample pairs of male subjects

library(tidyverse)

unrel_id <- read_tsv("data/1kgp/unrelated_subj.tsv") %>% pull(SAMPLE_NAME)

subj_info <- read_csv("data/1kgp/subject_info.csv") %>%
  filter(SAMPLE_NAME %in% unrel_id, sex == 1) #unrelated males

subj_info2 <- read_csv("data/1kgp/subject_info.csv") %>%
  filter(SAMPLE_NAME %in% unrel_id)


populations <- unique(subj_info$POPULATION)
eur_pops <- subj_info %>%
  filter(SUPER == "EUR") %>%
  pull(POPULATION) %>%
  unique()

afr_pops <- subj_info %>%
  filter(SUPER == "AFR") %>%
  pull(POPULATION) %>%
  unique()

# sample_pair <- function(df, populations){
#   s_pop <- sample(populations, 1)
#   subjects <- df %>%
#     filter(POPULATION == s_pop) %>%
#     pull(SAMPLE_NAME)
#   s_subj <- sample(subjects, 2)
#   return(c(s_pop, s_subj))
# }

sample_mixed_pair <- function(df, pop1, pop2){
  subjects1 <- df %>%
    filter(SUPER == pop1) %>%
    pull(SAMPLE_NAME)
  subjects2 <- df %>%
    filter(SUPER == pop2) %>%
    pull(SAMPLE_NAME)
  s_subj1 <- sample(subjects1, 1)
  s_subj2 <- sample(subjects2, 1)
  return(c(s_subj1, s_subj2))
}

set.seed(1071)
sampled_pairs <- replicate(200, sample_mixed_pair(subj_info, "AFR", "EUR")) %>% t() %>% as.data.frame()

write_csv(sampled_pairs, "data/sample_pairs_admx_8nov22.csv", col_names = FALSE, quote = "none")
