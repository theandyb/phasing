library(tidyverse)

populations <- read_tsv("data/1kgp/20131219.populations.tsv") %>%
  select(`Population Code`, `Super Population`) %>%
  rename(POPULATION = `Population Code`,
         SUPER = `Super Population`)

unrel <- read_tsv("data/1kgp/unrelated_subj.tsv") %>%
  select(SAMPLE_NAME, POPULATION) %>%
  left_join(populations, by = "POPULATION")


rel <- read_tsv("data/1kgp/related_subj.tsv") %>%
  select(SAMPLE_NAME, POPULATION) %>%
  left_join(populations, by = "POPULATION")

ped <- read_delim("data/1kgp/1kGP.3202_samples.pedigree_info.txt", delim = " ") %>%
  rename(SAMPLE_NAME = sampleID)

df <- bind_rows(unrel, rel) %>%
  full_join(ped, by = "SAMPLE_NAME")

write_csv(df, "data/1kgp/subject_info.csv", quote = "none")
