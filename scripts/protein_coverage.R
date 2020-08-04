library(tidyverse)
library(RSQLite)
library(unixtools)
library(data.table)
library(maditr)
library(lvplot)

source("libs/libs.R")

wdir <- "/vol/cloud/antonio/geogenetics/plant-arctic/results/mg_TO"

db <- "data/traits-search.sqlite"

con <- RSQLite::dbConnect(RSQLite::SQLite(), db)

# Calculate protein coverage ----------------------------------------------

# Viridiplantae -----------------------------------------------------------
busco_viridiplantae <- tbl(con, "busco_viridiplantae") %>%
  collect()

busco_viridiplantae_multi <- resolve_multi_mappings(busco_viridiplantae)

busco_viridiplantae_multi_covs <- busco_viridiplantae_multi %>%
  as.data.table() %>%
  dt_select(label, theader, tstart, tend, tlen) %>%
  setNames(c("label", "protein", "start", "end", "len")) %>%
  let(strand = "+") %>%
  arrange(protein, start, end) %>%
  setkeyv(., c("label", "protein"))

n_groups <- busco_viridiplantae_multi_covs %>% 
  dt_select(label, protein) %>% 
  unique() %>% 
  nrow()

pb <- txtProgressBar(min = 0, max = n_groups, style = 3)
busco_viridiplantae_multi_covs <- busco_viridiplantae_multi_covs[, {
  setTxtProgressBar(pb, .GRP);
  calculate_subject_coverage_dt(.SD);
}, by = list(label, protein), .SDcols = c("start", "end", "len")] 
close(pb)

# Fungi -------------------------------------------------------------------
busco_fungi <- tbl(con, "busco_fungi") %>%
  collect()

busco_fungi_multi <- resolve_multi_mappings(busco_fungi)

busco_fungi_multi_covs <- busco_fungi_multi %>%
  as.data.table() %>%
  dt_select(label, theader, tstart, tend, tlen) %>%
  setNames(c("label", "protein", "start", "end", "len")) %>%
  let(strand = "+") %>%
  arrange(protein, start, end) %>%
  setkeyv(., c("label", "protein"))

n_groups <- busco_fungi_multi_covs %>% 
  dt_select(label, protein) %>% 
  unique() %>% 
  nrow()

pb <- txtProgressBar(min = 0, max = n_groups, style = 3)
busco_fungi_multi_covs <- busco_fungi_multi_covs[, {
  setTxtProgressBar(pb, .GRP);
  calculate_subject_coverage_dt(.SD);
}, by = list(label, protein), .SDcols = c("start", "end", "len")] 
close(pb)

# Microbial traits --------------------------------------------------------
protein_traits <- tbl(con, "protein_traits") %>%
  collect()

protein_traits_multi <- resolve_multi_mappings(protein_traits)

protein_traits_multi_covs <- protein_traits_multi %>%
  as.data.table() %>%
  dt_select(label, theader, tstart, tend, tlen) %>%
  setNames(c("label", "protein", "start", "end", "len")) %>%
  let(strand = "+") %>%
  arrange(protein, start, end) %>%
  setkeyv(., c("label", "protein"))

n_groups <- protein_traits_multi_covs %>% 
  dt_select(label, protein) %>% 
  unique() %>% 
  nrow()

pb <- txtProgressBar(min = 0, max = n_groups, style = 3)
protein_traits_multi_covs <- protein_traits_multi_covs[, {
  setTxtProgressBar(pb, .GRP);
  calculate_subject_coverage_dt(.SD);
}, by = list(label, protein), .SDcols = c("start", "end", "len")] 
close(pb)


# Plotting ----------------------------------------------------------------
busco_viridiplantae_plots <- plot_alns(busco_viridiplantae_multi_covs, title = "Viriplantae")
busco_fungi_plots <- plot_alns(busco_fungi_multi_covs, title = "Fungi")
protein_traits_plots <- plot_alns(protein_traits_multi_covs, title = "Traits")

sample_stats <- tbl(con, "sample_stats") %>%
  select(label, num_seqs) %>%
  collect()

busco_viridiplantae_reads_mapped <- busco_viridiplantae_multi_covs %>%
  select(label, reads) %>%
  group_by(label) %>%
  count(sort = TRUE) %>%
  ungroup() %>%
  inner_join(sample_stats) %>%
  mutate(prop_reads_mapped = n/num_seqs,
         class = "Viridiplantae")

busco_fungi_reads_mapped <- busco_fungi_multi_covs %>%
  select(label, reads) %>%
  group_by(label) %>%
  count(sort = TRUE) %>%
  ungroup() %>%
  inner_join(sample_stats) %>%
  mutate(prop_reads_mapped = n/num_seqs,
         class = "Fungi")

protein_traits_reads_mapped <- protein_traits_multi_covs %>%
  select(label, reads) %>%
  group_by(label) %>%
  count(sort = TRUE) %>%
  ungroup() %>%
  inner_join(sample_stats) %>%
  mutate(prop_reads_mapped = n/num_seqs,
         class = "Traits")

busco_viridiplantae_reads_mapped %>%
  bind_rows(busco_fungi_reads_mapped) %>%
  bind_rows(protein_traits_reads_mapped) %>%
  mutate(class = fct_relevel(class, c("Traits", "Fungi", "Viridiplantae"))) %>%
  ggplot(aes(class, prop_reads_mapped)) +
  geom_lv(size = 0.5, width.method = "height", color = "#404040", width = 0.5, alpha = 1) +
  scale_y_log10(labels = scales::percent) +
  ggpubr::rotate() +
  theme_light() +
  xlab("Reads mapped") +
  ylab("")

# Saving tables -----------------------------------------------------------
# Search results
write_tsv(x = busco_viridiplantae, path = "results/busco-viridiplantae_search-results.tsv.gz")
write_tsv(x = busco_fungi, path = "results/busco-fungi_search-results.tsv.gz")
write_tsv(x = protein_traits, path = "results/protein-traits_search-results.tsv.gz")

# Resolved multi mappings
write_tsv(x = busco_viridiplantae_multi, path = "results/busco-viridiplantae_search-results-multi.tsv.gz")
write_tsv(x = busco_fungi_multi, path = "results/busco-fungi_search-results-multi.tsv.gz")
write_tsv(x = protein_traits_multi, path = "results/protein-traits_search-results-multi.tsv.gz")

# Coverage using resolved multi mapping results
write_tsv(x = busco_viridiplantae_multi_covs, path = "results/busco-viridiplantae_search-results-multi-covs.tsv.gz")
write_tsv(x = busco_fungi_multi_covs, path = "results/busco-fungi_search-results-multi-covs.tsv.gz")
write_tsv(x = protein_traits_multi_covs, path = "results/protein-traits_search-results-multi-covs.tsv.gz")

dbDisconnect(con)
