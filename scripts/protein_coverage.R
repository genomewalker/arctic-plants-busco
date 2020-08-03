library(tidyverse)
library(RSQLite)
library(unixtools)

Sys.setenv(RENV_PATHS_CACHE = "/vol/cloud/antonio/tmp/renv/cache")
unixtools::set.tempdir(path.expand("/vol/cloud/antonio/tmp"))

wdir <- "/vol/cloud/antonio/geogenetics/plant-arctic/results/mg_TO"

# Create DB from results --------------------------------------------------
busco_viridiplantae_files <- list.files(file.path(wdir, "busco_viridiplantae-homologs"), pattern = "tsv$", full.names = TRUE)
busco_fungi_files <- list.files(file.path(wdir, "busco_fungi-homologs"), pattern = "tsv$", full.names = TRUE)
protein_traits_files <- list.files(file.path(wdir, "traits-homologs"), pattern = "tsv$", full.names = TRUE)

col_names <- c("label", "query", "theader", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen")

# Calculate protein coverage ----------------------------------------------

# Viridiplantae -----------------------------------------------------------
busco_viridiplantae <- map_dfr(busco_viridiplantae_files, read_tsv, col_names = col_names) %>%
  mutate(strand = ifelse(qstart < qend, "+", "-"),
         qstart1 = ifelse(qstart < qend, qstart, qend),
         qend = ifelse(qstart > qend, qstart, qend),
         qstart = qstart1) %>%
  select(-qstart1)

busco_viridiplantae_multi <- resolve_multi_mappings(busco_viridiplantae)

busco_viridiplantae_multi_covs <- busco_viridiplantae_multi %>%
  select(label, theader, tstart, tend, tlen) %>%
  rename(protein = theader, start = tstart, end = tend, len = tlen) %>%
  mutate(strand = "+") %>%
  arrange(protein, start, end) %>%
  group_by(label, protein) %>%
  group_modify(~calculate_subject_coverage(.)) %>%
  ungroup()


# Fungi -------------------------------------------------------------------

busco_fungi <- map_dfr(busco_fungi_files, read_tsv, col_names = col_names) %>%
  mutate(strand = ifelse(qstart < qend, "+", "-"),
         qstart1 = ifelse(qstart < qend, qstart, qend),
         qend = ifelse(qstart > qend, qstart, qend),
         qstart = qstart1) %>%
  select(-qstart1)

busco_fungi_multi <- resolve_multi_mappings(busco_fungi)

busco_fungi_multi_covs <- busco_fungi_multi %>%
  select(label, theader, tstart, tend, tlen) %>%
  rename(protein = theader, start = tstart, end = tend, len = tlen) %>%
  mutate(strand = "+") %>%
  arrange(protein, start, end) %>%
  group_by(label, protein) %>%
  group_modify(~calculate_subject_coverage(.)) %>%
  ungroup()


# Microbial traits --------------------------------------------------------

protein_traits <- map_dfr(protein_traits_files, read_tsv, col_names = col_names) %>%
  mutate(strand = ifelse(qstart < qend, "+", "-"),
         qstart1 = ifelse(qstart < qend, qstart, qend),
         qend = ifelse(qstart > qend, qstart, qend),
         qstart = qstart1) %>%
  select(-qstart1)

protein_traits_multi <- resolve_multi_mappings(protein_traits)

protein_traits_multi_covs <- protein_traits_multi %>%
  select(label, theader, tstart, tend, tlen) %>%
  rename(protein = theader, start = tstart, end = tend, len = tlen) %>%
  mutate(strand = "+") %>%
  arrange(protein, start, end) %>%
  group_by(label, protein) %>%
  group_modify(~calculate_subject_coverage(.)) %>%
  ungroup()

# Plotting ----------------------------------------------------------------

busco_viridiplantae_plots <- plot_alns(busco_viridiplantae_multi_covs, title = "Viriplantae")
busco_fungi_plots <- plot_alns(busco_fungi_multi_covs, title = "Fungi")
protein_traits_plots <- plot_alns(protein_traits_multi_covs, title = "Traits")
