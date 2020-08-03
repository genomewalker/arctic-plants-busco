require(dplyr)

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

famli_in_r <- function(X, scale = 0.9, iters = 10){
  suppressWarnings({
    suppressMessages({
      if ("iter" %in% names(X)) {
        iter <- X$iter %>% unique()
      }else(
        iter <- 0
      )
      
      cat(paste0("Running Iter: ", iter + 1))
      w_bits <- X %>% 
        select(all_of(col_names)) %>%
        group_by(query, label) %>%
        mutate(W = bits/sum(bits),
               n_aln = n()) 
      
      tot_bits <- X %>% 
        select(label, theader, bits) %>%
        group_by(label, theader) %>%
        summarise(tot_bits = sum(bits),
                  n_aln_prot = n())
      
      l_bits <- w_bits %>% 
        filter(n_aln > 1) %>%
        inner_join(tot_bits) %>%
        mutate(L = W * tot_bits) %>%
        group_by(label, query) %>% 
        mutate(max_L = max(L),
               Keep_L = ifelse(L < scale * max_L, FALSE, TRUE)) %>%
        ungroup()
      l_bits_sel <- l_bits  %>%
        ungroup() %>%
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

plot_alns <- function(X, title = NULL){
  mode_val <- estimate_mode(X$depth_ratio)
  X_mode <- X %>% 
    filter(depth_ratio < mode_val) 
  
  X_int <- X %>%
    select(label, protein, depth_median, depth_ratio) %>%
    mutate(ints = cut(depth_ratio, breaks= seq(0, max(ceiling(depth_ratio)), by = 0.1), include.lowest = TRUE, right = FALSE)) 
  
  int <- X_int %>%
    group_by(ints) %>%
    summarise(mean = mean(depth_median),
              sd = sd(depth_median)) %>%
    ungroup() %>%
    mutate(n = row_number()) %>%
    filter(!is.na(sd)) 
  
  plot_int <- ggplot(int, aes(n, sd)) +
    geom_line() +
    theme_light() +
    xlab("Depth CV interval") +
    ylab("Median depth SD")
  
  int <- int[with(int,diff(sd)/diff(n)) %>%
               enframe(name = "step", value = "diff") %>% 
               filter(diff == 0) %>% 
               slice(which.min(step)) %>%
               head(1) %>% .$step,]
  
  cutoff <- X_int %>%
    filter(ints == int$ints) %>%
    slice(which.min(depth_ratio)) %>%
    head(1) %>%
    .$depth_ratio
  
  X_cut <- X %>% 
    filter(depth_ratio < cutoff) 
  
  plot_mode <- ggplot() +
    geom_point(data = X, aes(depth_ratio, depth_median), shape = 21, size = 2, fill = "grey", color = "black") +
    geom_point(data = X_cut, aes(depth_ratio, depth_median), fill = "#70B9E4", shape = 21, size = 2, color = "black") +
    scale_x_log10() +
    theme_light() +
    theme(panel.background = element_blank()) +
    xlab("Depth coefficient of variation") +
    ylab("Median depth")
  
  protein_order <- X_cut %>%
    group_by(protein) %>%
    mutate(median_depth = (median(depth_mean)),
           median_depth_log = log10(median(depth_mean))) %>%
    ungroup() %>%
    select(protein, median_depth, median_depth_log) %>%
    distinct() %>%
    arrange(desc(median_depth_log))
  
  X_dm <- X_cut %>%
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
    theme_light() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank()) +
    xlab("Proteins") +
    ylab("Mean coverage depth") +
    ggtitle(title)
  
  list(mode = mode_val, cutoff_val = cutoff,  X_mode = X_mode, X_cut = X_cut, plot_mode = plot_mode, plot_box = plot_box, n_prots = nrow(protein_order), plot_int = plot_int)
}
