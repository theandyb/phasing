# Sample pairs of male subjects

library(tidyverse)

unrel_id <- read_tsv("data/1kgp/unrelated_subj.tsv") %>% pull(SAMPLE_NAME)
# unrel_id <- read_tsv("/net/snowwhite/home/beckandy/research/phasing/data/gnomAD/1kgp_samples.txt", col_names = "SAMPLE_NAME") %>%
#   pull(SAMPLE_NAME)

subj_info <- read_csv("data/1kgp/subject_info.csv") %>%
  filter(SAMPLE_NAME %in% unrel_id, sex == 1) #unrelated males

#temp: generate male id list for vcf subsetting
# subj_info %>%
#   select(SAMPLE_NAME) %>%
#   write_tsv("/net/snowwhite/home/beckandy/research/phasing/data/male_ids.txt", col_names = FALSE)

# temp: generate sex file for fixploidy
# subj_info %>%
#   mutate(sex = "M") %>%
#   select(SAMPLE_NAME, sex) %>%
#   write_tsv("/net/snowwhite/home/beckandy/research/phasing/data/gnomAD/sex_list.txt", col_names = FALSE)

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

amr_pops <- subj_info %>%
  filter(SUPER == "AMR") %>%
  pull(POPULATION) %>%
  unique()

eas_pops <- subj_info %>%
  filter(SUPER == "EAS") %>%
  pull(POPULATION) %>%
  unique()

sas_pops <- subj_info %>%
  filter(SUPER == "SAS") %>%
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

set.seed(1990)
sampled_pairs_afr <- replicate(200, sample_pair(subj_info, afr_pops)) %>% t() %>% as.data.frame()
sampled_pairs_eur <- replicate(200, sample_pair(subj_info, eur_pops)) %>% t() %>% as.data.frame()
sampled_pairs_amr <- replicate(100, sample_pair(subj_info, amr_pops)) %>% t() %>% as.data.frame()
sampled_pairs_eas <- replicate(100, sample_pair(subj_info, eas_pops)) %>% t() %>% as.data.frame()
sampled_pairs_sas <- replicate(100, sample_pair(subj_info, sas_pops)) %>% t() %>% as.data.frame()

while(dim(sampled_pairs_afr %>% distinct())[1] < 200){
  n_sample <- 200 - dim(sampled_pairs_afr %>% distinct())[1]
  sampled_pairs_afr <- bind_rows(sampled_pairs_afr,
                                 replicate(n_sample, sample_pair(subj_info, afr_pops)) %>% t() %>% as.data.frame()) %>%
    distinct()
}

while(dim(sampled_pairs_eur %>% distinct())[1] < 200){
  n_sample <- 200 - dim(sampled_pairs_eur %>% distinct())[1]
  sampled_pairs_eur <- bind_rows(sampled_pairs_eur,
                                 replicate(n_sample, sample_pair(subj_info, eur_pops)) %>% t() %>% as.data.frame()) %>%
    distinct()
}

while(dim(sampled_pairs_amr %>% distinct())[1] < 100){
  n_sample <- 100 - dim(sampled_pairs_amr %>% distinct())[1]
  sampled_pairs_amr <- bind_rows(sampled_pairs_amr,
                                 replicate(n_sample, sample_pair(subj_info, amr_pops)) %>% t() %>% as.data.frame()) %>%
    distinct()
}

while(dim(sampled_pairs_eas %>% distinct())[1] < 100){
  n_sample <- 100 - dim(sampled_pairs_eas %>% distinct())[1]
  sampled_pairs_eas <- bind_rows(sampled_pairs_eas,
                                 replicate(n_sample, sample_pair(subj_info, eas_pops)) %>% t() %>% as.data.frame()) %>%
    distinct()
}

while(dim(sampled_pairs_sas %>% distinct())[1] < 100){
  n_sample <- 100 - dim(sampled_pairs_sas %>% distinct())[1]
  sampled_pairs_sas <- bind_rows(sampled_pairs_sas,
                                 replicate(n_sample, sample_pair(subj_info, sas_pops)) %>% t() %>% as.data.frame()) %>%
    distinct()
}

df <- bind_rows(sampled_pairs_eur, sampled_pairs_afr) %>%
  bind_rows(sampled_pairs_amr) %>%
  bind_rows(sampled_pairs_eas) %>%
  bind_rows(sampled_pairs_sas)
write_csv(df, "data/sample_pairs_9aug2023.csv", col_names = FALSE, quote = "none")
