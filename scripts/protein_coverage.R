library(tidyverse)
library(RSQLite)
library(data.table)
library(maditr)
library(lvplot)
library(gghalves)

source("libs/libs.R")

db <- "data/traits-search.sqlite"

con <- RSQLite::dbConnect(RSQLite::SQLite(), db)

# Define plotting theme
ggthemr::ggthemr(layout = "scientific", palette = "flat")

# Calculate protein coverage ----------------------------------------------

# Viridiplantae -----------------------------------------------------------
busco_viridiplantae <- tbl(con, "busco_viridiplantae") %>%
  collect() %>%
  filter(bits >= 60)

busco_viridiplantae_multi <- resolve_multi_mappings(busco_viridiplantae)

busco_viridiplantae_multi_covs <- busco_viridiplantae_multi %>%
  as.data.table() %>%
  dt_select(label, theader, tstart, tend, tlen) %>%
  setNames(c("label", "protein", "start", "end", "len")) %>%
  let(strand = "+") %>%
  setkeyv(., c("label", "protein"))

busco_viridiplantae_multi_covs <- get_protein_coverage(busco_viridiplantae_multi_covs)

# Fungi -------------------------------------------------------------------
busco_fungi <- tbl(con, "busco_fungi") %>%
  collect() %>%
  filter(bits >= 60)

busco_fungi_multi <- resolve_multi_mappings(busco_fungi)

busco_fungi_multi_covs <- busco_fungi_multi %>%
  as.data.table() %>%
  dt_select(label, theader, tstart, tend, tlen) %>%
  setNames(c("label", "protein", "start", "end", "len")) %>%
  let(strand = "+") %>%
  setkeyv(., c("label", "protein"))

busco_fungi_multi_covs <- get_protein_coverage(busco_fungi_multi_covs)

# Microbial traits --------------------------------------------------------
protein_traits <- tbl(con, "protein_traits") %>%
  collect() %>%
  filter(bits >= 60)

protein_traits_multi <- resolve_multi_mappings(protein_traits)

protein_traits_multi_covs <- protein_traits_multi %>%
  as.data.table() %>%
  dt_select(label, theader, tstart, tend, tlen) %>%
  setNames(c("label", "protein", "start", "end", "len")) %>%
  let(strand = "+") %>%
  setkeyv(., c("label", "protein"))

protein_traits_multi_covs <- get_protein_coverage(protein_traits_multi_covs)

# Plotting ----------------------------------------------------------------
busco_viridiplantae_plots <- process_alignments(busco_viridiplantae_multi_covs, title = "Viriplantae")
busco_fungi_plots <- process_alignments(busco_fungi_multi_covs, title = "Fungi")
protein_traits_plots <- process_alignments(protein_traits_multi_covs, title = "Traits")

# Basic stats -------------------------------------------------------------

types <- c("Traits", "Fungi", "Viridiplantae")
names(types) <- c("traits", "fungi", "viridiplantae")

types_colors <- c("#E85A5B", "#a3806d", "#a6e18f")
names(types_colors) <- c("Traits", "Fungi", "Viridiplantae")

# How many reads did we recruit

sample_stats <- tbl(con, "sample_stats") %>%
  select(label, num_seqs) %>%
  collect()


busco_viridiplantae_reads_mapped <- get_reads_stats(busco_viridiplantae_multi_covs, type = "Viridiplantae")
busco_fungi_reads_mapped <- get_reads_stats(busco_viridiplantae_multi_covs, type = "Fungi")
protein_traits_reads_mapped <- get_reads_stats(busco_viridiplantae_multi_covs, type = "Traits")

read_sample_plot <- busco_viridiplantae_reads_mapped %>%
  bind_rows(busco_fungi_reads_mapped) %>%
  bind_rows(protein_traits_reads_mapped) %>%
  mutate(type = fct_relevel(type, (types))) %>%
  ggplot(aes(type, prop_reads_mapped, fill = type)) +
  geom_half_boxplot(fill = "#1D1E22",
                    color = "black",
                    alpha = 0.4,
                    errorbar.draw = TRUE, 
                    width = 0.2, 
                    errorbar.length = 0.2, outlier.shape = NA) +
  geom_half_point_panel(transformation = ggbeeswarm:::position_quasirandom(width = 0.1), 
                        shape = 21, 
                        color = "black",
                        alpha = 0.9,
                        size = 2) +  
  scale_y_log10(labels = scales::percent) +
  scale_fill_manual(values = types_colors) +
  ggpubr::rotate() +
  theme(legend.position = "none",
        panel.background = element_blank()) +
  ylab("Reads mapped per sample") +
  xlab("")

# How many proteins have a certain coverage

total_proteins <- tbl(con, "proteins") %>%
  group_by(type) %>%
  count(name = "n_proteins") %>%
  collect() %>%
  mutate(type = recode(type, !!!types))


viridiplantae_prot_stats <- get_protein_stats(busco_viridiplantae_multi_covs, type = "Viridiplantae")
fungi_prot_stats <- get_protein_stats(busco_fungi_multi_covs, type = "Fungi")
traits_prot_stats <- get_protein_stats(protein_traits_multi_covs, type = "Traits")

protein_counts_plot <- viridiplantae_prot_stats$protein_counts %>%
  bind_rows(fungi_prot_stats$protein_counts) %>%
  bind_rows(traits_prot_stats$protein_counts) %>%
  inner_join(total_proteins) %>%
  mutate(prop_proteins = counts/n_proteins) %>%
  mutate(type = fct_relevel(type, types)) %>%
  ggplot(aes(type, prop_proteins, fill = type)) +
  geom_half_boxplot(fill = "#1D1E22",
                    color = "black",
                    alpha = 0.4,
                    errorbar.draw = TRUE, 
                    width = 0.2, 
                    errorbar.length = 0.2, outlier.shape = NA) +
  geom_half_point_panel(transformation = ggbeeswarm:::position_quasirandom(width = 0.1), 
                        shape = 21, 
                        color = "black",
                        alpha = 0.9,
                        size = 2) +  
  scale_y_log10(labels = scales::percent) +
  scale_fill_manual(values = types_colors) +
  ggpubr::rotate() +
  theme(legend.position = "none",
        panel.background = element_blank()) +
  ylab("Proteins found per sample") +
  xlab("")

# How many proteins are shared
colors <- c("#D0D5D9", "#CBC6C1", "#F8ECE0")
names(colors) <- c("Solid", "Shared", "Weak")
bind_rows(list(get_protein_sets(busco_viridiplantae_plots, type = "Viridiplantae"),
               get_protein_sets(busco_fungi_plots, type = "Fungi"),
               get_protein_sets(protein_traits_plots, type = "Traits"))) %>%
  pivot_longer(-type) %>%
  group_by(type) %>%
  mutate(total = sum(value),
         prop = value/total) %>%
  ungroup() %>%
  mutate(type = fct_relevel(type, types),
         name = fct_relevel(name, c("Solid", "Shared", "Weak"))) %>%
  ggplot(aes(type, prop, fill = name)) +
  geom_col(color = "black", size = 0.5) +
  ggpubr::rotate() +
  scale_fill_manual(values = colors, name = NULL, guide = guide_legend(reverse = TRUE)) +
  theme(legend.position = "top",
        line = element_blank(), rect = element_blank(),
        axis.ticks.x = element_line(),
        axis.text.y = element_text(margin = margin(r = -15)),
        legend.key.size = unit(1,"line")) +
  scale_y_continuous(labels = scales::percent) +
  xlab("") +
  ylab("Proportion of proteins")



# Saving plots ------------------------------------------------------------
ggpubr::ggarrange(busco_viridiplantae_plots$plot_depth_ratio_cutoff + ggtitle("Viridiplantae"),
                  busco_fungi_plots$plot_depth_ratio_cutoff + ggtitle("Fungi"),
                  protein_traits_plots$plot_depth_ratio_cutoff + ggtitle("Traits"),
                  nrow = 3, align = "hv")

ggpubr::ggarrange(busco_viridiplantae_plots$plot_intervals_depth + ggtitle("Viridiplantae"),
                  busco_fungi_plots$plot_intervals_depth + ggtitle("Fungi"),
                  protein_traits_plots$plot_intervals_depth + ggtitle("Traits"),
                  nrow = 3, align = "hv")

ggpubr::ggarrange(busco_viridiplantae_plots$plot_boxplots + ggtitle("Viridiplantae"),
                  busco_fungi_plots$plot_boxplots + ggtitle("Fungi"),
                  protein_traits_plots$plot_boxplots + ggtitle("Traits"),
                  nrow = 3, align = "hv")

ggpubr::ggarrange(read_sample_plot, protein_counts_plot, grid::nullGrob(), nrow = 1, align = "hv", widths = c(0.5, 0.5, 0.1))

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
write_tsv(x = busco_viridiplantae_multi_covs, path = "results/busco-viridiplantae_search-results-multi-covs_all.tsv.gz")
write_tsv(x = busco_fungi_multi_covs, path = "results/busco-fungi_search-results-multi-covs_all.tsv.gz")
write_tsv(x = protein_traits_multi_covs, path = "results/protein-traits_search-results-multi-covs_all.tsv.gz")

# Coverage using resolved multi mapping results - SOLID
write_tsv(x = busco_viridiplantae_plots$alns_filtered_depth_ratio_cutoff_solid, path = "results/busco-viridiplantae_search-results-multi-covs_solid.tsv.gz")
write_tsv(x = busco_fungi_plots$alns_filtered_depth_ratio_cutoff_solid, path = "results/busco-fungi_search-results-multi-covs_solid.tsv.gz")
write_tsv(x = protein_traits_plots$alns_filtered_depth_ratio_cutoff_solid, path = "results/protein-traits_search-results-multi-covs_solid.tsv.gz")

# Coverage using resolved multi mapping results - WEAK
write_tsv(x = busco_viridiplantae_plots$alns_filtered_depth_ratio_cutoff_weak, path = "results/busco-viridiplantae_search-results-multi-covs_weak.tsv.gz")
write_tsv(x = busco_fungi_plots$alns_filtered_depth_ratio_cutoff_weak, path = "results/busco-fungi_search-results-multi-covs_weak.tsv.gz")
write_tsv(x = protein_traits_plots$alns_filtered_depth_ratio_cutoff_weak, path = "results/protein-traits_search-results-multi-covs_weak.tsv.gz")

# Solid proteins with largest global median depths
write_tsv(x = busco_viridiplantae_plots$solid_protein_mediandepth_ranking, path = "results/busco-viridiplantae_top-mediandepth-proteins.tsv.gz")
write_tsv(x = busco_fungi_plots$solid_protein_mediandepth_ranking, path = "results/busco-fungi_top-mediandepth-proteins.tsv.gz")
write_tsv(x = protein_traits_plots$solid_protein_mediandepth_ranking, path = "results/protein-traits_top-mediandepth-proteins.tsv.gz")


dbDisconnect(con)
