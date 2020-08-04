library(RSQLite)
library(tidyverse)
library(janitor)

Sys.setenv(RENV_PATHS_CACHE = "/vol/cloud/antonio/tmp/renv/cache")
unixtools::set.tempdir(path.expand("/vol/cloud/antonio/tmp"))

wdir <- "/vol/cloud/antonio/geogenetics/plant-arctic/results/mg_TO"

db <- "data/traits-search.sqlite"

con <- RSQLite::dbConnect(RSQLite::SQLite(), db)

# Create DB from results --------------------------------------------------
busco_viridiplantae_files <- list.files(file.path(wdir, "busco_viridiplantae-homologs"), pattern = "tsv$", full.names = TRUE)
busco_fungi_files <- list.files(file.path(wdir, "busco_fungi-homologs"), pattern = "tsv$", full.names = TRUE)
protein_traits_files <- list.files(file.path(wdir, "traits-homologs"), pattern = "tsv$", full.names = TRUE)

col_names <- c("label", "query", "theader", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen")


# Get search results ------------------------------------------------------
# coverage:0.6 cov-mode:2 spaced-kmer-mode:0 sens:5.7 evalue:1e-20 sub-mat:PFAM30
# Viridiplantae -----------------------------------------------------------
busco_viridiplantae <- map_dfr(busco_viridiplantae_files, read_tsv, col_names = col_names) %>%
  mutate(strand = ifelse(qstart < qend, "+", "-"),
         qstart1 = ifelse(qstart < qend, qstart, qend),
         qend = ifelse(qstart > qend, qstart, qend),
         qstart = qstart1) %>%
  select(-qstart1)

copy_to(con, busco_viridiplantae, "busco_viridiplantae",
        temporary = FALSE,
        indexes = list(
          "label", "query", "theader"),
        overwrite = TRUE)

# Fungi -------------------------------------------------------------------
busco_fungi <- map_dfr(busco_fungi_files, read_tsv, col_names = col_names) %>%
  mutate(strand = ifelse(qstart < qend, "+", "-"),
         qstart1 = ifelse(qstart < qend, qstart, qend),
         qend = ifelse(qstart > qend, qstart, qend),
         qstart = qstart1) %>%
  select(-qstart1)

copy_to(con, busco_fungi, "busco_fungi",
        temporary = FALSE,
        indexes = list(
          "label", "query", "theader"),
        overwrite = TRUE)

# Microbial traits --------------------------------------------------------
protein_traits <- map_dfr(protein_traits_files, read_tsv, col_names = col_names) %>%
  mutate(strand = ifelse(qstart < qend, "+", "-"),
         qstart1 = ifelse(qstart < qend, qstart, qend),
         qend = ifelse(qstart > qend, qstart, qend),
         qstart = qstart1) %>%
  select(-qstart1)

copy_to(con, protein_traits, "protein_traits",
        temporary = FALSE,
        indexes = list(
          "label", "query", "theader"),
        overwrite = TRUE)

# List of samples
samples <- read_tsv("/vol/cloud/antonio/geogenetics/plant-arctic/mg_TO_wf/assets/samples-1M.txt")
copy_to(con, samples, "samples",
        temporary = FALSE,
        indexes = list(
          "label"),
        overwrite = TRUE)
# List of proteins
fungi_proteins <- read_tsv("/vol/cloud/antonio/geogenetics/plant-arctic/DB/TRAITS/fungi_ancestral_busco-db.headers", col_name = "protein") %>%
  mutate(type = "fungi")
viridiplantae_proteins <- read_tsv("/vol/cloud/antonio/geogenetics/plant-arctic/DB/TRAITS/viridiplantae_ancestral_busco-db.headers", col_name = "protein") %>%
  mutate(type = "viridiplantae")
traits_proteins <- read_tsv("/vol/cloud/antonio/geogenetics/plant-arctic/DB/TRAITS/traits_proteins-derep-db.headers", col_name = "protein") %>%
  mutate(type = "traits")

copy_to(con, fungi_proteins %>% bind_rows(viridiplantae_proteins) %>% bind_rows(traits_proteins), "proteins",
        temporary = FALSE,
        indexes = list(
          "protein"),
        overwrite = TRUE)

# Sample stats
stat_files <- list.files(file.path(wdir, "stats"), pattern = "txt$", full.names = TRUE)

stats <- map_dfr(stat_files, read_tsv, col_names = TRUE) %>%
 clean_names() %>%
  mutate(file = gsub("/vol/cloud/antonio/geogenetics/plant-arctic/mg_TO_wf/fastq/", "", file),
         file = gsub(".gz$", "", file)) %>%
  rename(label = file) %>%
  filter(label %in% samples$label)

copy_to(con, stats, "sample_stats",
        temporary = FALSE,
        indexes = list(
          "label"),
        overwrite = TRUE)

dbDisconnect(con)
