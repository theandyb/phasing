#COMMON FUNCTIONS

'%!in%' <- function(x,y)!('%in%'(x,y))

#q: why not just pull the start position of switches identified to have been the start of a flip?
get_flip_pos <- function(df){
  flip_list <- c()
  in_progress <- FALSE
  for(i in 1:(length(df$pos_start)-1)){
    if(df$pos_end[i] == df$pos_start[i+1]){
      if(!in_progress){
        in_progress <- TRUE
        flip_list <- c(flip_list, df$pos_end[i])
      } else{
        in_progress <- FALSE
      }
    } else{
      in_progress <- FALSE
    }
  }
  return(flip_list)
}

get_all_whatshap <- function(algo_name = "eagle",
                             whatshap_dir = "/net/snowwhite/home/beckandy/research/phasing/output/final_switch_errors/whatshap/",
                             n = 400,
                             pop = c(rep("EUR", 200), rep("AFR", 200))){
  results <- vector(mode = "list", length = n)
  for(i in 1:n){
    df <- read_tsv(paste0(whatshap_dir, "/", algo_name, "/eval_", i, ".tsv"), show_col_types = FALSE) %>%
      janitor::clean_names() %>%
      select(all_assessed_pairs, all_switches, all_switch_rate, all_switchflips, all_switchflip_rate) %>%
      mutate(switches = as.numeric(str_split(all_switchflips, "/")[[1]][1]),
             flips = as.numeric(str_split(all_switchflips, "/")[[1]][2])) %>%
      rename(n_het = all_assessed_pairs, total_errors = all_switches, error_rate = all_switch_rate) %>%
      select(-all_switchflips)
    df$id <- i
    results[[i]] <- df
  }
  final <- bind_rows(results)
  final$pop <- pop
  return(final)
}

switch_summary <- function(pair_id, eagle_dir, beagle_dir, shapeit_dir, gc_content_1kb, het_loc_dir, bin_size = 1000){
  switch_err_eagle <- read_csv(paste0(eagle_dir, "switch_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(bin_id = ceiling(pos_start / bin_size))
  switch_err_beagle <- read_csv(paste0(beagle_dir, "switch_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(bin_id = ceiling(pos_start / bin_size))
  switch_err_shapeit <- read_csv(paste0(shapeit_dir, "switch_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(bin_id = ceiling(pos_start / bin_size))

  het_df <- read_tsv(paste0(paste0(het_loc_dir, "pair_", pair_id, "_het_loc.txt")),
                     col_names = c("chr", "pos_start", "gt"), show_col_types = FALSE) %>%
    select(chr, pos_start) %>%
    mutate(pos_end = lead(pos_start)) %>%
    drop_na()

  switch_err_eagle <- switch_err_eagle %>%
    left_join({gc_content_1kb %>% select(bin_id, GC)}, by = "bin_id")

  switch_err_beagle <- switch_err_beagle %>%
    left_join({gc_content_1kb %>% select(bin_id, GC)}, by = "bin_id")

  switch_err_shapeit <- switch_err_shapeit %>%
    left_join({gc_content_1kb %>% select(bin_id, GC)}, by = "bin_id")

  # get positions of flips
  flip_pos_eagle <- get_flip_pos(switch_err_eagle)
  flip_pos_beagle <- get_flip_pos(switch_err_beagle)
  flip_pos_shapeit <- get_flip_pos(switch_err_shapeit)

  # Assign switches flip status
  switch_err_eagle$is_flip <- (switch_err_eagle$pos_start %in% flip_pos_eagle) |
    (switch_err_eagle$pos_end %in% flip_pos_eagle)

  switch_err_beagle$is_flip <- (switch_err_beagle$pos_start %in% flip_pos_beagle) |
    (switch_err_beagle$pos_end %in% flip_pos_beagle)

  switch_err_shapeit$is_flip <- (switch_err_shapeit$pos_start %in% flip_pos_shapeit) |
    (switch_err_shapeit$pos_end %in% flip_pos_shapeit)

  # Start of flip
  switch_err_eagle$start_flip <- (switch_err_eagle$pos_end %in% flip_pos_eagle)

  switch_err_beagle$start_flip <- (switch_err_beagle$pos_end %in% flip_pos_beagle)

  switch_err_shapeit$start_flip <- (switch_err_shapeit$pos_end %in% flip_pos_shapeit)

  # Append switch/flip locations to het_df
  het_df <- het_df %>%
    mutate(beagle_switch = pos_start %in% {
      switch_err_beagle %>%
        filter(!is_flip) %>%
        pull(pos_start)},
      beagle_flip = pos_start %in% {
        switch_err_beagle %>%
          filter(start_flip) %>%
          pull(pos_start)
      })
  het_df <- het_df %>%
    mutate(eagle_switch = pos_start %in% {
      switch_err_eagle %>%
        filter(!is_flip) %>%
        pull(pos_start)},
      eagle_flip = pos_start %in% {
        switch_err_eagle %>%
          filter(start_flip) %>%
          pull(pos_start)
      })
  het_df <- het_df %>%
    mutate(shapeit_switch = pos_start %in% {
      switch_err_shapeit %>%
        filter(!is_flip) %>%
        pull(pos_start)},
      shapeit_flip = pos_start %in% {
        switch_err_shapeit %>%
          filter(start_flip) %>%
          pull(pos_start)
      })
  het_df$id <- 1:length(het_df$chr)

  # stats we want to pull
  n_switch_eagle <- length(switch_err_eagle$pos_start)
  n_switch_beagle <- length(switch_err_beagle$pos_start)
  n_switch_shapeit <- length(switch_err_shapeit$pos_start)

  n_flip_eagle <- sum(switch_err_eagle$is_flip) / 2
  n_flip_beagle <- sum(switch_err_beagle$is_flip) / 2
  n_flip_shapeit <- sum(switch_err_shapeit$is_flip) / 2

  n_other_eagle <- n_switch_eagle - 2*n_flip_eagle
  n_other_beagle <- n_switch_beagle - 2*n_flip_beagle
  n_other_shapeit <- n_switch_shapeit - 2*n_flip_shapeit

  n_switch_cpg_eagle <- switch_err_eagle %>%
    filter(cpg_start == 1) %>%
    filter(is_flip != TRUE | start_flip == TRUE) %>%
    pull(is_flip) %>%
    length()
  n_switch_cpg_beagle <- switch_err_beagle %>%
    filter(cpg_start == 1) %>%
    filter(is_flip != TRUE | start_flip == TRUE) %>%
    pull(is_flip) %>%
    length()
  n_switch_cpg_shapeit <- switch_err_shapeit %>%
    filter(cpg_start == 1) %>%
    filter(is_flip != TRUE | start_flip == TRUE) %>%
    pull(is_flip) %>%
    length()

  n_flip_cpg_eagle <- switch_err_eagle %>%
    filter(cpg_start == 1) %>%
    pull(start_flip) %>%
    sum()
  n_flip_cpg_beagle <- switch_err_beagle %>%
    filter(cpg_start == 1) %>%
    pull(start_flip) %>%
    sum()
  n_flip_cpg_shapeit <- switch_err_shapeit %>%
    filter(cpg_start == 1) %>%
    pull(start_flip) %>%
    sum()

  n_other_cpg_eagle <- switch_err_eagle %>%
    filter(is_flip != TRUE) %>%
    pull(cpg_start) %>%
    sum()

  n_other_cpg_beagle <- switch_err_beagle %>%
    filter(is_flip != TRUE) %>%
    pull(cpg_start) %>%
    sum()

  n_other_cpg_shapeit <- switch_err_shapeit %>%
    filter(is_flip != TRUE) %>%
    pull(cpg_start) %>%
    sum()


  mean_gc_switch_eagle <- mean(switch_err_eagle$GC)
  mean_gc_switch_beagle <- mean(switch_err_beagle$GC)
  mean_gc_switch_shapeit <- mean(switch_err_shapeit$GC)

  # distance metrics
  median_dist_beagle <- switch_err_beagle %>%
    filter(is_flip == F | start_flip == T) %>% # filter out latter halves of flips
    pull(pos_start) %>%
    diff(lag = 1) %>%
    median()
  median_dist_eagle <- switch_err_eagle %>%
    filter(is_flip == F | start_flip == T) %>%
    pull(pos_start) %>%
    diff(lag = 1) %>%
    median()
  median_dist_shapeit <- switch_err_shapeit %>%
    filter(is_flip == F | start_flip == T) %>%
    pull(pos_start) %>%
    diff(lag = 1) %>%
    median()

  mean_dist_beagle <- switch_err_beagle %>%
    filter(is_flip == F | start_flip == T) %>% # filter out latter halves of flips
    pull(pos_start) %>%
    diff(lag = 1) %>%
    mean()
  mean_dist_eagle <- switch_err_eagle %>%
    filter(is_flip == F | start_flip == T) %>%
    pull(pos_start) %>%
    diff(lag = 1) %>%
    mean()
  mean_dist_shapeit <- switch_err_shapeit %>%
    filter(is_flip == F | start_flip == T) %>%
    pull(pos_start) %>%
    diff(lag = 1) %>%
    mean()

  # number hets between switches
  med_hets_switch_beagle <- het_df %>%
    filter(beagle_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  med_hets_switch_eagle <- het_df %>%
    filter(eagle_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  med_hets_switch_shapeit <- het_df %>%
    filter(shapeit_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  # number hets between flips
  med_hets_flip_beagle <- het_df %>%
    filter(beagle_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  med_hets_flip_eagle <- het_df %>%
    filter(eagle_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  med_hets_flip_shapeit <- het_df %>%
    filter(shapeit_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  # mean number of hets between switches and flips
  mean_hets_switch_beagle <- het_df %>%
    filter(beagle_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  mean_hets_switch_eagle <- het_df %>%
    filter(eagle_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  mean_hets_switch_shapeit <- het_df %>%
    filter(shapeit_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  # number hets between flips
  mean_hets_flip_beagle <- het_df %>%
    filter(beagle_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  mean_hets_flip_eagle <- het_df %>%
    filter(eagle_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  mean_hets_flip_shapeit <- het_df %>%
    filter(shapeit_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)



  return(data.frame(pair_id = pair_id,
                    n_switch_eagle = n_switch_eagle,
                    n_switch_beagle = n_switch_beagle,
                    n_switch_shapeit = n_switch_shapeit,
                    n_flip_eagle = n_flip_eagle,
                    n_flip_beagle = n_flip_beagle,
                    n_flip_shapeit = n_flip_shapeit,
                    n_other_eagle = n_other_eagle,
                    n_other_beagle = n_other_beagle,
                    n_other_shapeit = n_other_shapeit,
                    n_switch_cpg_eagle = n_switch_cpg_eagle,
                    n_switch_cpg_beagle = n_switch_cpg_beagle,
                    n_switch_cpg_shapeit = n_switch_cpg_shapeit,
                    n_flip_cpg_eagle = n_flip_cpg_eagle,
                    n_flip_cpg_beagle = n_flip_cpg_beagle,
                    n_flip_cpg_shapeit = n_flip_cpg_shapeit,
                    n_other_cpg_eagle = n_other_cpg_eagle,
                    n_other_cpg_beagle = n_other_cpg_beagle,
                    n_other_cpg_shapeit = n_other_cpg_shapeit,
                    mean_gc_switch_eagle = mean_gc_switch_eagle,
                    mean_gc_switch_beagle = mean_gc_switch_beagle,
                    mean_gc_switch_shapeit = mean_gc_switch_shapeit,
                    median_dist_shapeit = median_dist_shapeit,
                    median_dist_eagle = median_dist_eagle,
                    median_dist_beagle = median_dist_beagle,
                    mean_dist_eagle = mean_dist_eagle,
                    mean_dist_shapeit = mean_dist_shapeit,
                    mean_dist_beagle = mean_dist_beagle,
                    med_hets_switch_beagle = med_hets_switch_beagle,
                    med_hets_switch_eagle = med_hets_switch_eagle,
                    med_hets_switch_shapeit = med_hets_switch_shapeit,
                    med_hets_flip_beagle = med_hets_flip_beagle,
                    med_hets_flip_eagle = med_hets_flip_eagle,
                    med_hets_flip_shapeit = med_hets_flip_shapeit,
                    mean_hets_switch_beagle = mean_hets_switch_beagle,
                    mean_hets_switch_eagle = mean_hets_switch_eagle,
                    mean_hets_switch_shapeit = mean_hets_switch_shapeit,
                    mean_hets_flip_beagle = mean_hets_flip_beagle,
                    mean_hets_flip_eagle = mean_hets_flip_eagle,
                    mean_hets_flip_shapeit = mean_hets_flip_shapeit,
                    n_hets = length(het_df$pos_start)
  ))
}

switch_summary2 <- function(pair_id, eagle_dir, beagle_dir, shapeit_dir, gc_content_1kb, het_loc_dir, bin_size = 1000){
  switch_err_eagle <- read_csv(paste0(eagle_dir, "switch_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(bin_id = ceiling(pos_start / bin_size))
  switch_err_beagle <- read_csv(paste0(beagle_dir, "switch_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(bin_id = ceiling(pos_start / bin_size))
  switch_err_shapeit <- read_csv(paste0(shapeit_dir, "switch_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(bin_id = ceiling(pos_start / bin_size))

  # pos,gt,cpg,motif
  het_df <- read_csv(paste0(het_loc_dir, "pair_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(pos_end = lead(pos)) %>%
    rename(pos_start = pos) %>%
    drop_na()

  switch_err_eagle <- switch_err_eagle %>%
    left_join({gc_content_1kb %>% select(bin_id, GC)}, by = "bin_id")

  switch_err_beagle <- switch_err_beagle %>%
    left_join({gc_content_1kb %>% select(bin_id, GC)}, by = "bin_id")

  switch_err_shapeit <- switch_err_shapeit %>%
    left_join({gc_content_1kb %>% select(bin_id, GC)}, by = "bin_id")

  # get positions of flips
  flip_pos_eagle <- get_flip_pos(switch_err_eagle)
  flip_pos_beagle <- get_flip_pos(switch_err_beagle)
  flip_pos_shapeit <- get_flip_pos(switch_err_shapeit)

  # Assign switches flip status
  switch_err_eagle$is_flip <- (switch_err_eagle$pos_start %in% flip_pos_eagle) |
    (switch_err_eagle$pos_end %in% flip_pos_eagle)

  switch_err_beagle$is_flip <- (switch_err_beagle$pos_start %in% flip_pos_beagle) |
    (switch_err_beagle$pos_end %in% flip_pos_beagle)

  switch_err_shapeit$is_flip <- (switch_err_shapeit$pos_start %in% flip_pos_shapeit) |
    (switch_err_shapeit$pos_end %in% flip_pos_shapeit)

  # Start of flip
  switch_err_eagle$start_flip <- (switch_err_eagle$pos_end %in% flip_pos_eagle)

  switch_err_beagle$start_flip <- (switch_err_beagle$pos_end %in% flip_pos_beagle)

  switch_err_shapeit$start_flip <- (switch_err_shapeit$pos_end %in% flip_pos_shapeit)

  # Append switch/flip locations to het_df
  het_df <- het_df %>%
    mutate(beagle_switch = pos_start %in% {
      switch_err_beagle %>%
        filter(!is_flip) %>%
        pull(pos_start)},
      beagle_flip = pos_start %in% {
        switch_err_beagle %>%
          filter(start_flip) %>%
          pull(pos_start)
      })
  het_df <- het_df %>%
    mutate(eagle_switch = pos_start %in% {
      switch_err_eagle %>%
        filter(!is_flip) %>%
        pull(pos_start)},
      eagle_flip = pos_start %in% {
        switch_err_eagle %>%
          filter(start_flip) %>%
          pull(pos_start)
      })
  het_df <- het_df %>%
    mutate(shapeit_switch = pos_start %in% {
      switch_err_shapeit %>%
        filter(!is_flip) %>%
        pull(pos_start)},
      shapeit_flip = pos_start %in% {
        switch_err_shapeit %>%
          filter(start_flip) %>%
          pull(pos_start)
      })
  het_df$id <- 1:length(het_df$pos_start)

  # stats we want to pull
  n_switch_eagle <- length(switch_err_eagle$pos_start)
  n_switch_beagle <- length(switch_err_beagle$pos_start)
  n_switch_shapeit <- length(switch_err_shapeit$pos_start)

  n_flip_eagle <- sum(switch_err_eagle$is_flip) / 2
  n_flip_beagle <- sum(switch_err_beagle$is_flip) / 2
  n_flip_shapeit <- sum(switch_err_shapeit$is_flip) / 2

  n_other_eagle <- n_switch_eagle - 2*n_flip_eagle
  n_other_beagle <- n_switch_beagle - 2*n_flip_beagle
  n_other_shapeit <- n_switch_shapeit - 2*n_flip_shapeit

  n_switch_cpg_eagle <- switch_err_eagle %>%
    filter(cpg_start == 1) %>%
    filter(is_flip != TRUE | start_flip == TRUE) %>%
    pull(is_flip) %>%
    length()
  n_switch_cpg_beagle <- switch_err_beagle %>%
    filter(cpg_start == 1) %>%
    filter(is_flip != TRUE | start_flip == TRUE) %>%
    pull(is_flip) %>%
    length()
  n_switch_cpg_shapeit <- switch_err_shapeit %>%
    filter(cpg_start == 1) %>%
    filter(is_flip != TRUE | start_flip == TRUE) %>%
    pull(is_flip) %>%
    length()

  n_flip_cpg_eagle <- switch_err_eagle %>%
    filter(cpg_start == 1) %>%
    pull(start_flip) %>%
    sum()
  n_flip_cpg_beagle <- switch_err_beagle %>%
    filter(cpg_start == 1) %>%
    pull(start_flip) %>%
    sum()
  n_flip_cpg_shapeit <- switch_err_shapeit %>%
    filter(cpg_start == 1) %>%
    pull(start_flip) %>%
    sum()

  n_other_cpg_eagle <- switch_err_eagle %>%
    filter(is_flip != TRUE) %>%
    pull(cpg_start) %>%
    sum()

  n_other_cpg_beagle <- switch_err_beagle %>%
    filter(is_flip != TRUE) %>%
    pull(cpg_start) %>%
    sum()

  n_other_cpg_shapeit <- switch_err_shapeit %>%
    filter(is_flip != TRUE) %>%
    pull(cpg_start) %>%
    sum()


  mean_gc_switch_eagle <- mean(switch_err_eagle$GC)
  mean_gc_switch_beagle <- mean(switch_err_beagle$GC)
  mean_gc_switch_shapeit <- mean(switch_err_shapeit$GC)

  # distance metrics
  median_dist_beagle <- switch_err_beagle %>%
    filter(is_flip == F | start_flip == T) %>% # filter out latter halves of flips
    pull(pos_start) %>%
    diff(lag = 1) %>%
    median()
  median_dist_eagle <- switch_err_eagle %>%
    filter(is_flip == F | start_flip == T) %>%
    pull(pos_start) %>%
    diff(lag = 1) %>%
    median()
  median_dist_shapeit <- switch_err_shapeit %>%
    filter(is_flip == F | start_flip == T) %>%
    pull(pos_start) %>%
    diff(lag = 1) %>%
    median()

  mean_dist_beagle <- switch_err_beagle %>%
    filter(is_flip == F | start_flip == T) %>% # filter out latter halves of flips
    pull(pos_start) %>%
    diff(lag = 1) %>%
    mean()
  mean_dist_eagle <- switch_err_eagle %>%
    filter(is_flip == F | start_flip == T) %>%
    pull(pos_start) %>%
    diff(lag = 1) %>%
    mean()
  mean_dist_shapeit <- switch_err_shapeit %>%
    filter(is_flip == F | start_flip == T) %>%
    pull(pos_start) %>%
    diff(lag = 1) %>%
    mean()

  # number hets between switches
  med_hets_switch_beagle <- het_df %>%
    filter(beagle_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  med_hets_switch_eagle <- het_df %>%
    filter(eagle_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  med_hets_switch_shapeit <- het_df %>%
    filter(shapeit_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  # number hets between flips
  med_hets_flip_beagle <- het_df %>%
    filter(beagle_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  med_hets_flip_eagle <- het_df %>%
    filter(eagle_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  med_hets_flip_shapeit <- het_df %>%
    filter(shapeit_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    median(na.rm = T)

  # mean number of hets between switches and flips
  mean_hets_switch_beagle <- het_df %>%
    filter(beagle_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  mean_hets_switch_eagle <- het_df %>%
    filter(eagle_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  mean_hets_switch_shapeit <- het_df %>%
    filter(shapeit_switch) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  # number hets between flips
  mean_hets_flip_beagle <- het_df %>%
    filter(beagle_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  mean_hets_flip_eagle <- het_df %>%
    filter(eagle_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  mean_hets_flip_shapeit <- het_df %>%
    filter(shapeit_flip) %>%
    mutate(next_id = lead(id)) %>%
    mutate(hets_btw = next_id - id - 1) %>%
    pull(hets_btw) %>%
    mean(na.rm = T)

  prop_het_cpg <- sum(het_df$cpg) / length(het_df$cpg)



  return(data.frame(pair_id = pair_id,
                    n_switch_eagle = n_switch_eagle,
                    n_switch_beagle = n_switch_beagle,
                    n_switch_shapeit = n_switch_shapeit,
                    n_flip_eagle = n_flip_eagle,
                    n_flip_beagle = n_flip_beagle,
                    n_flip_shapeit = n_flip_shapeit,
                    n_other_eagle = n_other_eagle,
                    n_other_beagle = n_other_beagle,
                    n_other_shapeit = n_other_shapeit,
                    n_switch_cpg_eagle = n_switch_cpg_eagle,
                    n_switch_cpg_beagle = n_switch_cpg_beagle,
                    n_switch_cpg_shapeit = n_switch_cpg_shapeit,
                    n_flip_cpg_eagle = n_flip_cpg_eagle,
                    n_flip_cpg_beagle = n_flip_cpg_beagle,
                    n_flip_cpg_shapeit = n_flip_cpg_shapeit,
                    n_other_cpg_eagle = n_other_cpg_eagle,
                    n_other_cpg_beagle = n_other_cpg_beagle,
                    n_other_cpg_shapeit = n_other_cpg_shapeit,
                    mean_gc_switch_eagle = mean_gc_switch_eagle,
                    mean_gc_switch_beagle = mean_gc_switch_beagle,
                    mean_gc_switch_shapeit = mean_gc_switch_shapeit,
                    median_dist_shapeit = median_dist_shapeit,
                    median_dist_eagle = median_dist_eagle,
                    median_dist_beagle = median_dist_beagle,
                    mean_dist_eagle = mean_dist_eagle,
                    mean_dist_shapeit = mean_dist_shapeit,
                    mean_dist_beagle = mean_dist_beagle,
                    med_hets_switch_beagle = med_hets_switch_beagle,
                    med_hets_switch_eagle = med_hets_switch_eagle,
                    med_hets_switch_shapeit = med_hets_switch_shapeit,
                    med_hets_flip_beagle = med_hets_flip_beagle,
                    med_hets_flip_eagle = med_hets_flip_eagle,
                    med_hets_flip_shapeit = med_hets_flip_shapeit,
                    mean_hets_switch_beagle = mean_hets_switch_beagle,
                    mean_hets_switch_eagle = mean_hets_switch_eagle,
                    mean_hets_switch_shapeit = mean_hets_switch_shapeit,
                    mean_hets_flip_beagle = mean_hets_flip_beagle,
                    mean_hets_flip_eagle = mean_hets_flip_eagle,
                    mean_hets_flip_shapeit = mean_hets_flip_shapeit,
                    n_hets = length(het_df$pos_start),
                    prop_het_cpg = prop_het_cpg
  ))
}

load_het_pos <- function(id, het_loc_dir){
  fname <- paste0(het_loc_dir, "pair_", id, ".csv")
  return(read_csv(fname, show_col_types = FALSE))
}

load_switch_df <- function(id, switch_dir){
  fname <- paste0(switch_dir, "switch_", id, ".csv")
  results <- read_csv(fname, show_col_types = FALSE)
  flip_loc <- get_flip_pos(results)
  results$is_flip <- (results$pos_start %in% flip_loc) | (results$pos_end %in% flip_loc)
  results$flip <- results$pos_end %in% flip_loc
  # results <- results %>%
  #   filter(flip | (!is_flip & !flip))
  # return(results %>% select(-is_flip))
  return(results)
}

annotate_het_list <- function(id, het_loc_dir, beagle_dir, eagle_dir, shapeit_dir, df_freq){
  result <- load_het_pos(id, het_loc_dir)
  df_b <- load_switch_df(id, beagle_dir)
  df_e <- load_switch_df(id, eagle_dir)
  df_s <- load_switch_df(id, shapeit_dir)

  result$switch_b <- result$pos %in% {df_b %>% filter(!is_flip) %>% pull(pos_start)}
  result$flip_b <- result$pos %in% {df_b %>% filter(flip) %>% pull(pos_start)}
  result$flip2_b <- result$pos %in% {df_b %>% filter(is_flip & !flip) %>% pull(pos_start)}

  result$switch_e <- result$pos %in% {df_e %>% filter(!is_flip) %>% pull(pos_start)}
  result$flip_e <- result$pos %in% {df_e %>% filter(flip) %>% pull(pos_start)}
  result$flip2_e <- result$pos %in% {df_e %>% filter(is_flip & !flip) %>% pull(pos_start)}

  result$switch_s <- result$pos %in% {df_s %>% filter(!is_flip) %>% pull(pos_start)}
  result$flip_s <- result$pos %in% {df_s %>% filter(flip) %>% pull(pos_start)}
  result$flip2_s <- result$pos %in% {df_s %>% filter(is_flip & !flip) %>% pull(pos_start)}

  result <- left_join(result,
                      {df_freq %>%
                          select(POS, maf)},
                      by = c("pos" = "POS"))

  return(result)
}

get_error_counts <- function(ids, het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, df_freq){
  n <- length(ids)
  results <- vector(mode = "list", length = 2)

  df <- annotate_het_list(ids[1], het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, df_freq)
  df_switch <- df %>%
    filter(switch_b | switch_e | switch_s) %>%
    mutate(beagle = ifelse(switch_b, "Y", "N"),
           eagle = ifelse(switch_e, "Y", "N"),
           shapeit = ifelse(switch_s, "Y", "N")) %>%
    mutate(count_cat = paste(beagle, eagle, shapeit, sep="_")) %>%
    pull(count_cat) %>%
    table() %>%
    as.data.frame() %>%
    rename(count_cat = ".")

  df_flip <- df %>%
    filter(flip_b | flip_e | flip_s) %>%
    mutate(beagle = ifelse(flip_b, "Y", "N"),
           eagle = ifelse(flip_e, "Y", "N"),
           shapeit = ifelse(flip_s, "Y", "N")) %>%
    mutate(count_cat = paste(beagle, eagle, shapeit, sep="_")) %>%
    pull(count_cat) %>%
    table() %>%
    as.data.frame() %>%
    rename(count_cat = ".")

  df_switch_noflip <- df %>%
    filter((switch_b | switch_e | switch_s) & !(flip_b | flip2_b | flip_e | flip2_e | flip_s | flip2_s)) %>%
    mutate(beagle = ifelse(switch_b, "Y", "N"),
           eagle = ifelse(switch_e, "Y", "N"),
           shapeit = ifelse(switch_s, "Y", "N")) %>%
    mutate(count_cat = paste(beagle, eagle, shapeit, sep="_")) %>%
    pull(count_cat) %>%
    table() %>%
    as.data.frame() %>%
    rename(count_cat = ".")

  for(i in 2:n){
    df <- annotate_het_list(ids[i], het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, df_freq)
    df_switch2 <- df %>%
      filter(switch_b | switch_e | switch_s) %>%
      mutate(beagle = ifelse(switch_b, "Y", "N"),
             eagle = ifelse(switch_e, "Y", "N"),
             shapeit = ifelse(switch_s, "Y", "N")) %>%
      mutate(count_cat = paste(beagle, eagle, shapeit, sep="_")) %>%
      pull(count_cat) %>%
      table() %>%
      as.data.frame() %>%
      rename(count_cat = ".", Freq2 = Freq)

    df_switch <- df_switch %>%
      full_join(df_switch2, by="count_cat") %>%
      replace_na(list("Freq" = 0, "Freq2" = 0)) %>%
      mutate(Freq = Freq + Freq2) %>%
      select(-Freq2)


    df_flip2 <- df %>%
      filter(flip_b | flip_e | flip_s) %>%
      mutate(beagle = ifelse(flip_b, "Y", "N"),
             eagle = ifelse(flip_e, "Y", "N"),
             shapeit = ifelse(flip_s, "Y", "N")) %>%
      mutate(count_cat = paste(beagle, eagle, shapeit, sep="_")) %>%
      pull(count_cat) %>%
      table() %>%
      as.data.frame() %>%
      rename(count_cat = ".", Freq2 = Freq)

    df_flip <- df_flip %>%
      full_join(df_flip2, by="count_cat") %>%
      replace_na(list("Freq" = 0, "Freq2" = 0)) %>%
      mutate(Freq = Freq + Freq2) %>%
      select(-Freq2)

    df_switch_noflip2 <- df %>%
      filter((switch_b | switch_e | switch_s) & !(flip_b | flip2_b | flip_e | flip2_e | flip_s | flip2_s)) %>%
      mutate(beagle = ifelse(switch_b, "Y", "N"),
             eagle = ifelse(switch_e, "Y", "N"),
             shapeit = ifelse(switch_s, "Y", "N")) %>%
      mutate(count_cat = paste(beagle, eagle, shapeit, sep="_")) %>%
      pull(count_cat) %>%
      table() %>%
      as.data.frame() %>%
      rename(count_cat = ".", Freq2 = Freq)

    df_switch_noflip <- df_switch_noflip %>%
      full_join(df_switch_noflip2, by="count_cat") %>%
      replace_na(list("Freq" = 0, "Freq2" = 0)) %>%
      mutate(Freq = Freq + Freq2) %>%
      select(-Freq2)
  }
  results[[1]] <- df_switch
  results[[2]] <- df_flip
  results[[3]] <- df_switch_noflip
  return(results)
}

count_to_ven <- function(x){
  df <- x %>%
    mutate(beagle = ifelse(str_sub(count_cat, 1, 1)=="Y", 1, 0),
           eagle = ifelse(str_sub(count_cat, 3, 3)=="Y", 1, 0),
           shapeit = ifelse(str_sub(count_cat, 5, 5)== "Y", 1, 0)) %>%
    select(beagle, eagle, shapeit, Freq)

  n_b <- df %>% filter(beagle == 1, eagle == 0, shapeit == 0) %>% pull(Freq)
  n_e <- df %>% filter(beagle == 0, eagle == 1, shapeit == 0) %>% pull(Freq)
  n_s <- df %>% filter(beagle == 0, eagle == 0, shapeit == 1) %>% pull(Freq)
  n_b_e <- df %>% filter(beagle == 1, eagle == 1, shapeit == 0) %>% pull(Freq)
  n_b_s <- df %>% filter(beagle == 1, eagle == 0, shapeit == 1) %>% pull(Freq)
  n_e_s <- df %>% filter(beagle == 0, eagle == 1, shapeit == 1) %>% pull(Freq)
  n_b_e_s <- df %>% filter(beagle == 1, eagle == 1, shapeit == 1) %>% pull(Freq)

  data.frame(beagle = c(rep(TRUE, n_b), rep(FALSE, n_e),rep(FALSE, n_s), rep(TRUE, n_b_e), rep(TRUE, n_b_s), rep(FALSE, n_e_s) , rep(TRUE, n_b_e_s)),
             eagle = c(rep(FALSE, n_b), rep(TRUE, n_e),rep(FALSE, n_s), rep(TRUE, n_b_e), rep(FALSE, n_b_s), rep(TRUE, n_e_s) , rep(TRUE, n_b_e_s)),
             shapeit = c(rep(FALSE, n_b), rep(FALSE, n_e), rep(TRUE, n_s), rep(FALSE, n_b_e), rep(TRUE, n_b_s), rep(TRUE, n_e_s) , rep(TRUE, n_b_e_s))
  ) %>%
    vennCounts() %>%
    vennDiagram() %>% print()
}

get_errors_mb <- function(ids, bin_size, het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, df_freq){

  final_df <- annotate_het_list(ids[1], het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, df_freq) %>%
    mutate(bin_id = ceiling(pos / bin_size)) %>%
    group_by(bin_id) %>%
    summarize(beagle_s = sum(switch_b),
              eagle_s = sum(switch_e),
              shapeit_s = sum(switch_s),
              beagle_f = sum(flip_b),
              eagle_f = sum(flip_e),
              shapeit_f = sum(flip_s),
              hets = n())

  for(i in 2:length(ids)){
    df2 <- annotate_het_list(ids[i], het_loc_dir, beagle_switch_dir, eagle_switch_dir, shapeit_switch_dir, df_freq) %>%
      mutate(bin_id = ceiling(pos / bin_size)) %>%
      group_by(bin_id) %>%
      summarize(beagle_s2 = sum(switch_b),
                eagle_s2 = sum(switch_e),
                shapeit_s2 = sum(switch_s),
                beagle_f2 = sum(flip_b),
                eagle_f2 = sum(flip_e),
                shapeit_f2 = sum(flip_s),
                hets2 = n())
    final_df <- full_join(final_df, df2, by="bin_id") %>%
      replace_na(list("beagle_s" = 0, "eagle_s" = 0, "shapeit_s" = 0,
                      "beagle_s2" = 0, "eagle_s2" = 0, "shapeit_s2" = 0,
                      "beagle_f" = 0, "eagle_f" = 0, "shapeit_f" = 0,
                      "beagle_f2" = 0, "eagle_f2" = 0, "shapeit_f2" = 0,
                      "hets" = 0, "hets2" = 0)) %>%
      mutate(beagle_s = beagle_s + beagle_s2,
             eagle_s = eagle_s + eagle_s2,
             shapeit_s = shapeit_s + shapeit_s2,
             beagle_f = beagle_f + beagle_f2,
             eagle_f = eagle_f + eagle_f2,
             shapeit_f = shapeit_f + shapeit_f2,
             hets = hets + hets2) %>%
      select(bin_id, beagle_s, eagle_s, shapeit_s, beagle_f, eagle_f, shapeit_f, hets)
  }
  return(final_df %>% arrange(bin_id))
}

