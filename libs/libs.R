library(dplyr)
library(data.table)
library(maditr)
# Modification from https://www.biorxiv.org/content/10.1101/295352v2
create_aln_vector <- function(start, end, len) {
  v <- vector(length = len)
  v[start:end] <- 1
  v
}

calculate_subject_coverage <- function(X){
  X <- X %>% select(start, end, len)
  subject_len <- X$len %>% unique() 
  m <- pmap(X, create_aln_vector) %>% 
    simplify2array() %>% t()
  depths <- colSums(m) 
  
  tibble(cov = length(depths[depths != 0])/subject_len,
         depth_mean = mean(depths),
         depth_median = median(depths),
         depth_sd = sd(depths),
         depth_ratio = depth_sd/depth_mean,
         reads = nrow(X))
}

calculate_subject_coverage_dt <- function(X){
  depths = pmap(X, create_aln_vector) %>% simplify2array() %>% t() %>% colSums();
  subject_len = X$len %>% unique();
  cov = length(depths[depths != 0])/subject_len;
  depth_mean = mean(depths);
  depth_median = median(depths);
  depth_sd = sd(depths);
  depth_ratio = depth_sd/depth_mean;
  reads = nrow(X);
  list(cov = cov, depth_mean = depth_mean, depth_median = depth_median, depth_sd = depth_sd, depth_ratio = depth_ratio, reads = reads)
}

famli_in_r <- function(X, scale = 0.9, iters = 10){
  suppressWarnings({
    suppressMessages({
      if ("iter" %in% names(X)) {
        iter <- X$iter %>% unique()
      }else(
        iter <- 0
      )
      
      X <- as.data.table(X)
      
      key_cols <- c("label","query")
      setkeyv(X, key_cols)
      
      col_names <- c("label", "query", "theader", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen")
      
      cat(paste0("Running Iter: ", iter + 1))
      
      w_bits <- X %>% 
        dt_select(col_names) %>%
        let(W = bits/sum(bits),
            n_aln = .N, by = list(query, label))
      
      tot_bits <- X %>% 
        dt_select(label, theader, bits) %>%
        dt_summarise(tot_bits = sum(bits),
                     n_aln_prot = .N, by = list(label, theader))
      
      l_bits <- w_bits %>% 
        dt_filter(n_aln > 1) %>%
        dt_inner_join(tot_bits) %>%
        let(L = W * tot_bits) %>%
        let(max_L = max(L),
            Keep_L = ifelse(L < scale * max_L, FALSE, TRUE), by = list(label, query)) %>% as_tibble()
      
      l_bits_sel <- l_bits  %>%
        arrange(label, query) %>%
        filter(Keep_L == TRUE)
      
      keep_processing_n <- l_bits_sel %>% 
        group_by(label, query) %>% 
        add_count(name = "refined_n") %>% 
        filter(refined_n > 1) %>% nrow()
      
      cat(paste0(" Left: ", keep_processing_n, "... Filtering done\n") )
      
      keep_processing <- keep_processing_n == 0
      
      iter <- iter + 1
      
      if (!keep_processing & iter < iters){
        famli_in_r(X = l_bits_sel %>% mutate(iter = iter))
      }else{
        l_bits_sel
      }
    })
  })
}

resolve_multi_mappings <- function(alns, scale = 0.9, iters = 10){
  # Resolved multi mapping using ML approach from FAMLI
  suppressWarnings({
    suppressMessages({
      col_names <- c("label", "query", "theader", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen")
      multi_mapped <- famli_in_r(X = alns, scale = scale, iters = iters)
      alns %>%
        anti_join(multi_mapped %>% select(label, query) %>% distinct()) %>%
        bind_rows(multi_mapped %>% select(all_of(col_names)))
    })
  })
}


estimate_mode <- function(x) {
  if (length(x) == 2) {
    mean(x)
  } else if (length(x) == 1){
    x
  }else if (length(x) == 0){
    NA
  }else {
    d <- density(x)
    d$x[which.max(d$y)]
  }
}

process_alignments <- function(X, title = NULL){
  
  mode_val <- estimate_mode(X$depth_ratio)
  
  X_mode_solid <- X %>% 
    filter(depth_ratio < mode_val) 
  
  X_mode_weak <- X %>% 
    filter(depth_ratio >= mode_val) 
  
  X_int <- X %>%
    select(label, protein, depth_median, depth_ratio) %>%
    mutate(ints = cut(depth_ratio, breaks = seq(0, max(ceiling(depth_ratio)), by = 0.1), include.lowest = TRUE, right = FALSE)) 
  
  int <- X_int %>%
    group_by(ints) %>%
    summarise(mean = mean(depth_median),
              sd = sd(depth_median)) %>%
    ungroup() %>%
    mutate(n = row_number(),
           sd = ifelse(is.na(sd), 0, sd))
  
  int_diff <- with(int,diff(sd)/diff(n)) %>%
    enframe(name = "step", value = "diff")   
  
  int_zeroes <- int_diff %>% 
    filter(diff == 0) 
  
  find_longest_zeroes <- function(X, int_diff) {
    tibble(step = X, 
           median = median(int_diff %>% filter(step >= X) %>% .$diff),
           mean = mean(int_diff %>% filter(step >= X) %>% .$diff),
           n = length(int_diff %>% filter(step >= X) %>% .$diff))
  }
  
  idx <- map_dfr(int_zeroes$step, find_longest_zeroes, int_diff = int_diff) %>%
    filter(median == 0) %>%
    arrange(step) %>%
    head(1)
  
  if (is_empty(idx$step)){
    int_sel <- int[with(int,diff(sd)/diff(n)) %>%
                     enframe(name = "step", value = "diff") %>%
                     filter(diff == 0) %>%
                     slice(which.min(step)) %>%
                     head(1) %>% .$step,]
  }else{
    int_sel <- int[idx$step, ]
  }
  
  cutoff <- X_int %>%
    filter(ints == int_sel$ints) %>%
    slice(which.min(depth_ratio)) %>%
    head(1) %>%
    .$depth_ratio
  
  X_cut_solid <- X %>% 
    filter(depth_ratio < cutoff) 
  
  X_cut_weak <- X %>% 
    filter(depth_ratio >= cutoff) 
  
  plot_int <- ggplot(int, aes(n, sd)) +
    geom_vline(xintercept = idx$step, color = "#E85A5B", size = 1) +
    geom_line(color = "#454345", size = 1) +
    theme(panel.background = element_blank()) +
    xlab("Depth CV interval") +
    ylab("Median depth SD")
  
  plot_cut <- ggplot() +
    geom_point(data = X, aes(depth_ratio, depth_median), shape = 21, size = 2, fill = "grey", color = "black") +
    geom_point(data = X_cut_solid, aes(depth_ratio, depth_median), fill = "#70B9E4", shape = 21, size = 2, color = "black") +
    scale_x_log10() +
    theme(panel.background = element_blank()) +
    xlab("Depth coefficient of variation") +
    ylab("Median depth")
  
  protein_order <- X_cut_solid %>%
    group_by(protein) %>%
    mutate(median_depth = (median(depth_mean)),
           median_depth_log = log10(median(depth_mean))) %>%
    ungroup() %>%
    select(protein, median_depth, median_depth_log) %>%
    distinct() %>%
    arrange(desc(median_depth_log))
  
  X_dm <- X_cut_solid %>%
    mutate(protein = fct_relevel(protein, protein_order$protein)) %>%
    droplevels()
  
  formatBack <- function(x) scales::comma(10^x)   
  plot_box <- ggplot() +
    geom_boxplot(data = X_dm, aes(x = reorder(protein, -log10(depth_mean), FUN=median, order = TRUE), y = log10(depth_mean)), size = 0.2, color = "black", fill = "grey70", outlier.shape = 21, outlier.colour = "black", outlier.fill = "grey70", outlier.size = 0.8) +
    geom_hline(yintercept = log10(1), color = "#E85A5B") +
    stat_summary(data = X_dm, aes(x = reorder(protein, -log10(depth_mean), FUN=median, order = TRUE), y = log10(depth_mean)),
                 fun = median, 
                 geom = "line",
                 group = 1,
                 color = "#F5A428", 
                 size = 1) +
    #scale_y_continuous(labels = scales::trans_format("log10")) +
    scale_y_continuous(labels= formatBack) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank()) +
    xlab("Proteins") +
    ylab("Mean coverage depth") +
    ggtitle(title)
  
  list(depth_ratio_mode = mode_val, 
       depth_ratio_cutoff = cutoff,  
       alns_filtered_depth_ratio_mode_solid = X_mode_solid,
       alns_filtered_depth_ratio_mode_weak = X_mode_weak,
       alns_filtered_depth_ratio_cutoff_solid = X_cut_solid, 
       alns_filtered_depth_ratio_cutoff_weak = X_cut_weak, 
       plot_depth_ratio_cutoff = plot_cut, 
       plot_boxplots = plot_box, 
       n_proteins = nrow(protein_order), 
       plot_intervals_depth = plot_int, 
       depth_intervals = int, 
       depth_intervals_diff = int_diff,
       solid_protein_mediandepth_ranking = protein_order)
}

get_protein_coverage <- function(X){
  n_groups <- X %>% 
    dt_select(label, protein) %>% 
    unique() %>% 
    nrow()
  
  pb <- txtProgressBar(min = 0, max = n_groups, style = 3)
  X <- X[, {
    setTxtProgressBar(pb, .GRP);
    calculate_subject_coverage_dt(.SD);
  }, by = list(label, protein), .SDcols = c("start", "end", "len")] 
  close(pb)  
  return(X)
}

get_reads_stats <- function(X, type){
  X %>%
    select(label, reads) %>%
    group_by(label) %>%
    count(sort = TRUE) %>%
    ungroup() %>%
    inner_join(sample_stats) %>%
    mutate(prop_reads_mapped = n/num_seqs,
           type = type)
}

get_protein_stats <- function(X, type){
  protein_counts <- X %>%
    select(label, protein) %>%
    group_by(label) %>%
    count(name = "counts") %>%
    mutate(type = type) %>%
    ungroup()
  unique_proteins <- X$protein %>% uniqueN()
  list(protein_counts = protein_counts, unique_proteins = unique_proteins)
}

get_protein_sets <- function(X, type){
  solid <- X$alns_filtered_depth_ratio_cutoff_soli$protein %>% unique()
  weak <- X$alns_filtered_depth_ratio_cutoff_weak$protein %>% unique()
  tibble(
    type = type,
    Solid = setdiff(solid, weak) %>% length(),
    Weak = setdiff(weak, solid) %>% length(),
    Shared = intersect(solid, weak) %>% length())
}

