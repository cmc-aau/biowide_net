# Analysis of the networks

## Setting

### Loading libraries
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(ggrepel))
suppressMessages(library(wesanderson))
suppressMessages(library(pheatmap))

### Paths and settings
project.dir <- "/home/fdelogu/biowide_net"
data.dir <- file.path(project.dir, "data")
results.dir <- file.path(project.dir, "results")
figure.dir <- file.path(results.dir, "figures")
table.dir <- file.path(results.dir, "tables")
dir.create(figure.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table.dir, recursive = TRUE, showWarnings = FALSE)

# If an indicator index is not literally contained in each matching subnetwork
# name, replace the automatic matching below with an explicit condition lookup.
target.subnetwork <- "LateDryRich"

### Helper functions
read.data <- function(name, dir = data.dir) fread(file.path(dir, name)) %>% as_tibble()
save.plot <- function(plot, name, width, height, ...) {
  ggsave(file.path(figure.dir, name), plot = plot, dpi = 300,
         width = width, height = height, units = "cm", ...)
}
standardize <- function(x) {
  if (sum(!is.na(x)) < 2 || is.na(sd(x, na.rm = TRUE)) || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}
format.node <- function(x) {
  x %>% str_replace_all("sp\\.", "SP") %>% str_replace_all("\\.\\.", ";") %>%
    str_replace_all("\\.", " ") %>% str_replace_all(";", ". ") %>%
    str_replace_all("SP", "sp.")
}

### Loading data
prok.taxonomy <- read.data("2025-09-19_bg_tax_filt.tsv")
euk.taxonomy <- read.data("2025-09-19_ag_tax_filt.tsv")
prok.indicators <- read.data("2025-09-19_indicator_bg_rel.csv")
euk.indicators <- read.data("2025-09-19_indicator_ag_pa.csv")
node.metrics <- read.data("node_metrics.csv", file.path(results.dir, "topology"))
edge.metrics <- read.data("edge_metrics.csv", file.path(results.dir, "topology"))
net.metrics <- read.data("net_metrics.csv", file.path(results.dir, "topology"))

# Removed combined.tax: it was loaded but never used in the original analysis.

### Taxonomy and indicator lookup
# The original network IDs are retained for every join. Display formatting happens only at export.
tax.shortcut <- bind_rows(
  euk.taxonomy %>% transmute(Unit = str_replace_all(Species, " ", "."), Kingdom, Genus, Species),
  prok.taxonomy %>% transmute(Unit = genome, Kingdom = str_remove(Kingdom, "^d__"), Genus, Species)
) %>% distinct()

if (anyDuplicated(tax.shortcut$Unit)) {
  stop("Taxonomy contains repeated Unit IDs. Resolve these before joining to node/edge metrics.")
}

ind.shortcut <- bind_rows(
  prok.indicators %>% transmute(Unit = genome, index = as.character(index)),
  euk.indicators %>% transmute(Unit = str_replace_all(Species, " ", "."), index = as.character(index))
) %>% filter(!is.na(Unit), !is.na(index), index != "") %>% distinct()

# One row per (species, subnetwork). A match on ANY indicator condition gives status 1.
# This reproduces the original substring rule, but without rowwise joins or duplicated species.
subnetworks <- bind_rows(node.metrics %>% select(subnetwork),
                         edge.metrics %>% select(subnetwork),
                         net.metrics %>% select(subnetwork)) %>% distinct()
indicator.lookup <- crossing(Unit = unique(ind.shortcut$Unit), subnetwork = subnetworks$subnetwork) %>%
  left_join(ind.shortcut, by = "Unit", relationship = "many-to-many") %>%
  mutate(matched = map2_lgl(index, subnetwork, ~ !is.na(.x) && str_detect(.y, fixed(.x)))) %>%
  group_by(Unit, subnetwork) %>%
  summarise(indicator = as.integer(any(matched)), .groups = "drop")

# Independent lookup validation: exactly one status per Unit/subnetwork.
stopifnot(!anyDuplicated(indicator.lookup[c("Unit", "subnetwork")]))

### Shared node and edge metrics
node.analysis <- node.metrics %>%
  group_by(subnetwork) %>%
  mutate(positive_degree.std = standardize(positive_degree),
         positive_strength.std = standardize(positive_strength),
         betweenness.std = standardize(betweenness)) %>%
  ungroup() %>%
  left_join(tax.shortcut, by = c("node" = "Unit"), relationship = "many-to-one") %>%
  left_join(indicator.lookup, by = c("node" = "Unit", "subnetwork" = "subnetwork"), relationship = "many-to-one") %>%
  mutate(indicator = replace_na(indicator, 0L),
         `Indicator status` = if_else(indicator == 1, "Indicator", "Not indicator"),
         Species = if_else(Kingdom %in% c("Archaea", "Bacteria"), str_remove(Species, "^s__"), str_replace_all(Species, "\\.", " ")),
         Genus = if_else(Kingdom %in% c("Archaea", "Bacteria"), str_remove(Genus, "^g__"), Genus),
         taxonomy = case_when(
           Kingdom %in% c("Animalia", "Fungi", "Plantae") ~ Species,
           Kingdom %in% c("Archaea", "Bacteria") & !str_detect(coalesce(Species, ""), "novel") ~ Species,
           Kingdom %in% c("Archaea", "Bacteria") & str_detect(coalesce(Species, ""), "novel") ~ str_c(Genus, Species, sep = " ")))

ind.weight.tax <- edge.metrics %>%
  left_join(indicator.lookup %>% rename(source_node = Unit, source.indicator = indicator),
            by = c("source_node", "subnetwork"), relationship = "many-to-one") %>%
  left_join(indicator.lookup %>% rename(target_node = Unit, target.indicator = indicator),
            by = c("target_node", "subnetwork"), relationship = "many-to-one") %>%
  left_join(tax.shortcut %>% select(Unit, Kingdom) %>% rename(source_node = Unit, Source_kingdom = Kingdom),
            by = "source_node", relationship = "many-to-one") %>%
  left_join(tax.shortcut %>% select(Unit, Kingdom) %>% rename(target_node = Unit, Target_kingdom = Kingdom),
            by = "target_node", relationship = "many-to-one") %>%
  mutate(source.indicator = replace_na(source.indicator, 0L),
         target.indicator = replace_na(target.indicator, 0L),
         count.indicators = source.indicator + target.indicator,
         source_community = case_when(Source_kingdom %in% c("Animalia", "Fungi", "Plantae") ~ "above",
                                      Source_kingdom %in% c("Archaea", "Bacteria") ~ "below"),
         target_community = case_when(Target_kingdom %in% c("Animalia", "Fungi", "Plantae") ~ "above",
                                      Target_kingdom %in% c("Archaea", "Bacteria") ~ "below"),
         community = case_when(is.na(source_community) | is.na(target_community) ~ NA_character_,
                               source_community == target_community ~ "Intra", TRUE ~ "Between"))

stopifnot(nrow(node.analysis) == nrow(node.metrics), nrow(ind.weight.tax) == nrow(edge.metrics))
if (anyNA(node.analysis$Kingdom)) warning("Some nodes lack taxonomy matches.")
if (anyNA(ind.weight.tax$community)) warning("Some edges lack a community assignment.")

kingdom.color <- c("Animalia" = "#E41A1C", "Archaea" = "#FF7F00", "Bacteria" = "#006699",
                   "Fungi" = "#A6761D", "Plantae" = "#4DAF4A")
indicator.color <- c("Indicator" = "#000000", "Not indicator" = "#CCCCCC")
plot.colors <- c("Not indicator" = wes_palette("GrandBudapest1")[1],
                 "Indicator" = wes_palette("GrandBudapest2")[2])

## Exploratory

### PCA of network metrics
pca_metrics <- c("density", "positive_edge_fraction", "mean_positive_strength",
                 "mean_negative_strength", "positive_transitivity", "positive_modularity")
invalid.metrics <- pca_metrics[map_lgl(net.metrics[pca_metrics], ~ anyNA(.x) || !is.numeric(.x) || sd(.x) == 0)]
if (length(invalid.metrics)) message("Excluding PCA metrics containing NA/non-numeric/zero variance: ", paste(invalid.metrics, collapse = ", "))
pca_metrics <- setdiff(pca_metrics, invalid.metrics)
if (length(pca_metrics) < 2 || nrow(net.metrics) < 3) stop("PCA requires at least two variable metrics and three subnetworks.")
pca <- prcomp(net.metrics[pca_metrics], center = TRUE, scale. = TRUE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_scores <- tibble(subnetwork = net.metrics$subnetwork, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
loadings <- as_tibble(pca$rotation[, 1:2, drop = FALSE], rownames = "metric")
arrow_scale <- 0.8 * min(diff(range(pca_scores$PC1)) / diff(range(loadings$PC1)),
                          diff(range(pca_scores$PC2)) / diff(range(loadings$PC2)))
loadings <- loadings %>% mutate(across(c(PC1, PC2), ~ .x * arrow_scale))

pc12.loadings.plot <- ggplot(pca_scores, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = subnetwork)) +
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               inherit.aes = FALSE, color = "grey60", arrow = arrow(length = grid::unit(0.18, "cm"))) +
  geom_text_repel(data = loadings, aes(PC1, PC2, label = metric), inherit.aes = FALSE,
                  color = "grey60", size = 3.5) +
  labs(x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2])) +
  theme_classic(base_size = 13)
pc12.loadings.plot
save.plot(pc12.loadings.plot, "PC1_2_loadings_net_metrics.png", 15, 15)

mat <- net.metrics %>% select(all_of(pca_metrics)) %>% as.matrix()
rownames(mat) <- net.metrics$subnetwork
# A heatmap of standardized network metrics, not PCA loadings.
png(file.path(figure.dir, "PCA_loadings_heatmap.png"), width = 15, height = 15, units = "cm", res = 300)
pheatmap(mat, scale = "column", cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA)
dev.off()

### Node degree, strength and betweenness
deg.strength.plot <- node.analysis %>%
  ggplot(aes(positive_degree.std, positive_strength.std, color = subnetwork)) +
  geom_hline(yintercept = 2.5, color = "grey50", linetype = "dotted") +
  geom_vline(xintercept = 2.5, color = "grey50", linetype = "dotted") +
  geom_point() +
  geom_text_repel(data = ~ filter(.x, positive_degree.std > 2.5, positive_strength.std > 2.5),
                  aes(label = Genus), max.overlaps = 100) +
  theme_bw() + theme(aspect.ratio = 1)
deg.strength.plot
save.plot(deg.strength.plot, "deg_str.png", 25, 20)

plot.betweenness <- function(data, degree.cutoff = 2, betweenness.cutoff = 4, size = 2.5, overlaps = 35) {
  ggplot(data, aes(positive_degree.std, betweenness.std)) +
    geom_point(aes(fill = Kingdom, color = `Indicator status`), shape = 21, size = 2.5) +
    geom_hline(yintercept = betweenness.cutoff, color = "grey50", linetype = "dotted") +
    geom_vline(xintercept = degree.cutoff, color = "grey50", linetype = "dotted") +
    geom_text_repel(data = ~ filter(.x, positive_degree.std >= degree.cutoff,
                                   betweenness.std >= betweenness.cutoff),
                    aes(label = taxonomy), max.overlaps = overlaps, force = 100,
                    size = size, segment.color = "grey70", segment.alpha = 0.7) +
    scale_fill_manual(values = kingdom.color) +
    scale_color_manual(values = indicator.color) +
    theme_bw() + theme(aspect.ratio = 1) +
    labs(x = "Standardized positive degree", y = "Standardized node betweenness")
}

deg.betweenness.plot <- plot.betweenness(node.analysis) + facet_wrap(~subnetwork, nrow = 2)
deg.betweenness.plot
save.plot(deg.betweenness.plot, "deg_bet.png", 40, 20)
save.plot(deg.betweenness.plot, "deg_bet.svg", 40, 20)

deg.betweenness.LDR.plot <- node.analysis %>% filter(subnetwork == target.subnetwork) %>%
  plot.betweenness(degree.cutoff = 2, betweenness.cutoff = 4, size = 3, overlaps = 30)
deg.betweenness.LDR.plot
save.plot(deg.betweenness.LDR.plot, "deg_bet_LDR.png", 20, 20)
save.plot(deg.betweenness.LDR.plot, "deg_bet_LDR.svg", 20, 20)

### Positive edge weights across communities and indicator pairs
pd <- position_dodge(width = 0.8)
strength.indicators.communities.plot <- ind.weight.tax %>%
  filter(weight > 0, !is.na(community)) %>%
  mutate(community = recode(community, Between = "Between communities", Intra = "Intra communities"),
         count.indicators = factor(count.indicators, levels = 0:2,
                                   labels = c("No indicator species in pair", "One indicator species in pair",
                                              "Both indicator species in pair"))) %>%
  rename(Community = community) %>%
  ggplot(aes(subnetwork, weight, fill = Community)) +
  geom_violin(alpha = 0.5, position = pd) +
  geom_boxplot(width = 0.125, alpha = 0.5, position = pd) +
  stat_compare_means(method = "t.test", label = "p.signif", label.x.npc = 0, label.y.npc = 0.95,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf),
                                        symbols = c("***", "**", "*", "ns"))) +
  facet_wrap(~count.indicators, nrow = 1) +
  scale_fill_manual(values = c("Between communities" = wes_palette("GrandBudapest1")[1],
                               "Intra communities" = wes_palette("GrandBudapest2")[2])) +
  coord_flip(ylim = c(0, 1)) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1.1)) +
  labs(x = "Ecospace", y = "Modified Jaccard distance (J>0)")
strength.indicators.communities.plot
save.plot(strength.indicators.communities.plot, "strenght_ind_agbg.png", 20, 20)

### Isolate LateDryRich network and indicators
LateDryRich.net <- ind.weight.tax %>% filter(subnetwork == target.subnetwork) %>%
  mutate(across(c(source_node, target_node), format.node))

# Preserve one node per original network node; add a display label only in exported tables.
LateDryRich.nodes <- node.analysis %>% filter(subnetwork == target.subnetwork) %>%
  semi_join(tibble(node = c(edge.metrics$source_node[edge.metrics$subnetwork == target.subnetwork],
                              edge.metrics$target_node[edge.metrics$subnetwork == target.subnetwork])) %>% distinct(), by = "node") %>%
  mutate(node = format.node(node), Species = format.node(Species),
         Species = str_remove(Species, "^s__"),
         indicator = if_else(indicator == 1, "Yes", "No"))

fwrite(LateDryRich.net, file.path(table.dir, "LateDryRich_net.csv"))
fwrite(LateDryRich.nodes, file.path(table.dir, "LateDryRich_nodes.csv"))

### Analysis on LateDryRich network
LateDryRich.plot <- node.analysis %>%
  filter(subnetwork == target.subnetwork, Kingdom != "Archaea") %>%
  select(Kingdom, `Indicator status`, positive_degree, negative_degree,
         betweenness, positive_clustering) %>%
  pivot_longer(cols = -c(Kingdom, `Indicator status`), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric, betweenness = "Betweenness", negative_degree = "Degree (J<0)",
                         positive_degree = "Degree (J>0)", positive_clustering = "Clustering (J>0)"),
         Metric = factor(Metric, levels = c("Degree (J>0)", "Degree (J<0)", "Betweenness", "Clustering (J>0)")),
         `Indicator status` = factor(`Indicator status`, levels = c("Not indicator", "Indicator"))) %>%
  ggplot(aes(`Indicator status`, Value, color = `Indicator status`)) +
  geom_jitter(alpha = 0.25) + geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.125, alpha = 0.5) +
  stat_compare_means(method = "t.test", label = "p.signif", label.x.npc = 0.5, label.y.npc = 0.8,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf),
                                        symbols = c("***", "**", "*", "ns"))) +
  facet_wrap(Metric ~ Kingdom, scales = "free_y", nrow = 4) +
  scale_color_manual(values = plot.colors) +
  theme_bw() +
  theme(legend.position = "bottom", text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1.1)) +
  labs(x = "", y = "")
LateDryRich.plot
save.plot(LateDryRich.plot, "LateDryRich_subplot.png", 25, 25)
