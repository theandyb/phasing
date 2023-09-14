#COMMON FUNCTIONS

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

switch_summary <- function(pair_id, eagle_dir, beagle_dir, shapeit_dir, gc_content_1kb, bin_size = 1000){
  switch_err_eagle <- read_csv(paste0(eagle_dir, "switch_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(bin_id = ceiling(pos_start / bin_size))
  switch_err_beagle <- read_csv(paste0(beagle_dir, "switch_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(bin_id = ceiling(pos_start / bin_size))
  switch_err_shapeit <- read_csv(paste0(shapeit_dir, "switch_", pair_id, ".csv"), show_col_types = FALSE) %>%
    mutate(bin_id = ceiling(pos_start / bin_size))

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
                    mean_dist_beagle = mean_dist_beagle
  ))
}
