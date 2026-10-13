### Loading libraries

suppressMessages(library(tidyverse))

### Settings

base.dir <- "/home/fdelogu/biowide_net"
data.dir <- file.path(base.dir, "data")
topology.dir <- file.path(base.dir, "results", "topology")

### Loading node metrics

node.metrics <- read_csv(file.path(topology.dir, "node_metrics.csv"),
                         show_col_types = FALSE)

### Loading taxonomy

prok.taxonomy <- read_tsv(file.path(data.dir, "2025-09-19_bg_tax_filt.tsv"),
                          show_col_types = FALSE)
euk.taxonomy <- read_tsv(file.path(data.dir, "2025-09-19_ag_tax_filt.tsv"),
                         show_col_types = FALSE)

tax.shortcut <- bind_rows(
  euk.taxonomy %>%
    mutate(node = str_replace_all(Species, fixed(" "), ".")),
  prok.taxonomy %>%
    mutate(node = genome)
) %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species),
                ~ str_remove(.x, "^[a-z]__"))) %>%
  select(node, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  distinct(node, .keep_all = TRUE)

### Loading indicators

prok.indicators <- read_csv(file.path(data.dir, "2025-09-19_indicator_bg_rel.csv"),
                            show_col_types = FALSE)
euk.indicators <- read_csv(file.path(data.dir, "2025-09-19_indicator_ag_pa.csv"),
                           show_col_types = FALSE)

ind.shortcut <- bind_rows(
  prok.indicators %>%
    transmute(node = genome, index),
  euk.indicators %>%
    transmute(node = str_replace_all(Species, fixed(" "), "."), index)
) %>%
  distinct(node, index) %>%
  group_by(node) %>%
  summarise(n.index = n_distinct(index),
            indices = list(sort(unique(index))),
            .groups = "drop")

### Subnetwork conditions

network.conditions <- tibble(subnetwork = unique(node.metrics$subnetwork)) %>%
  mutate(succession = str_extract(subnetwork, "^(Early|Late)"),
         moisture = str_extract(subnetwork, "(Dry|Wet)"),
         fertility = str_extract(subnetwork, "(Poor|Rich)$")) %>%
  pivot_longer(cols = c(succession, moisture, fertility),
               names_to = "axis", values_to = "condition") %>%
  group_by(subnetwork) %>%
  summarise(conditions = list(condition), .groups = "drop")

stopifnot(!anyNA(unlist(network.conditions$conditions)))

### Loading OTU matrices and computing prevalence

subnetworks <- sort(unique(node.metrics$subnetwork))

prevalence.table <- map_dfr(subnetworks, function(net) {
  
  otu <- read.csv(file.path(base.dir, paste0("results/PA_split_OTUs/jaccard/", net, "_OTU.csv")),
                  row.names = 1, check.names = FALSE)
  
  tibble(subnetwork = net,
         node = str_replace_all(rownames(otu), fixed(" "), "."),
         prevalence = rowSums(as.matrix(otu) > 0),
         n_samples = ncol(otu)) %>%
    mutate(prevalence_frac = prevalence / n_samples)
  
})

### Building the annotated node metrics

node.metrics <- node.metrics %>%
  left_join(tax.shortcut, by = "node") %>%
  left_join(ind.shortcut, by = "node") %>%
  left_join(network.conditions, by = "subnetwork") %>%
  left_join(prevalence.table, by = c("subnetwork", "node")) %>%
  mutate(n.index = replace_na(n.index, 0L),
         indices = map(indices, ~ if (is.null(.x)) character(0) else .x),
         n.match = map2_int(indices, conditions, ~ sum(.x %in% .y)),
         matches = map2_chr(indices, conditions,
                            ~ paste(intersect(.x, .y), collapse = ";")),
         class = case_when(n.match > 0 ~ "Matching",
                           n.index > 0 ~ "Nonmatching",
                           TRUE ~ "Never indicator")) %>%
  select(-indices, -conditions)

### Quality checks

stopifnot(!anyDuplicated(node.metrics[c("subnetwork", "node")]))
stopifnot(nrow(node.metrics) == 18055)

node.metrics %>%
  count(subnetwork, class) %>%
  pivot_wider(names_from = class, values_from = n, values_fill = 0) %>%
  print(n = Inf)

node.metrics %>%
  summarise(missing.taxonomy = sum(is.na(Kingdom)),
            missing.prevalence = sum(is.na(prevalence_frac)),
            missing.phylum = sum(is.na(Phylum))) %>%
  print()

### Saving annotated table

write_csv(node.metrics,
          file.path(data.dir, "indicator_taxonomy_prevalence_audit.csv"))
