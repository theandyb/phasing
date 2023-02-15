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

sample_pair <- function(df, populations){
  s_pop <- sample(populations, 1)
  subjects <- df %>%
    filter(POPULATION == s_pop) %>%
    pull(SAMPLE_NAME)
  s_subj <- sample(subjects, 2)
  return(c(s_pop, s_subj))
}

set.seed(1071)
sampled_pairs_afr <- replicate(200, sample_pair(subj_info, afr_pops)) %>% t() %>% as.data.frame()
sampled_pairs_eur <- replicate(200, sample_pair(subj_info, eur_pops)) %>% t() %>% as.data.frame()

df <- bind_rows(sampled_pairs_eur, sampled_pairs_afr)
write_csv(df, "data/sample_pairs_16aug2022.csv", col_names = FALSE, quote = "none")
