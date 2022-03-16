# Sample pairs of male subjects

library(tidyverse)

unrel_id <- read_tsv("data/unrelated_subj.tsv") %>% pull(SAMPLE_NAME)

subj_info <- read_csv("data/subject_info.csv") %>%
  filter(SAMPLE_NAME %in% unrel_id, sex == 1) #unrelated males

populations <- unique(subj_info$POPULATION)
eur_pops <- subj_info %>%
  filter(SUPER == "EUR") %>%
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
sampled_pairs <- replicate(100, sample_pair(subj_info, eur_pops)) %>% t() %>% as.data.frame()
write_csv(sampled_pairs, "data/sample_pairs.csv", col_names = FALSE, quote = "none")
