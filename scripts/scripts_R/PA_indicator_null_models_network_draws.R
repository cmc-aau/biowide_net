# Indicator taxa and network topology
# Complete-profile, taxon-aware null model; networks remain fixed.
# Requires: tidyverse; future and future.apply only for parallel permutations.

### Settings
suppressMessages(library(tidyverse))

base.dir <- "/home/fdelogu/biowide_net"
data.dir <- file.path(base.dir, "data")   # Folder containing the three input CSV files
out.dir <- file.path(base.dir, "results", "indicator_null_models")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

n.perm <- 1999L      # Start with 199; increase to 9999 for the final analysis
n.cores <- 2L        # Set to 1 to disable parallel processing
chunk.size <- 100L  # Permutations per parallel task
seed <- 77081L
specifications <- c(3L, 5L) # Prevalence tertiles (primary) and quintiles (sensitivity)

### Loading and checking data
node.metrics <- readr::read_csv(file.path(data.dir, "indicator_taxonomy_prevalence_audit.csv"), show_col_types = FALSE)
euk.indicators <- readr::read_csv(file.path(data.dir, "2025-09-19_indicator_ag_pa.csv"), show_col_types = FALSE)
prok.indicators <- readr::read_csv(file.path(data.dir, "2025-09-19_indicator_bg_rel.csv"), show_col_types = FALSE)

metrics <- c("positive_degree", "positive_strength", "negative_degree", "betweenness")
required <- c("node", "subnetwork", "Kingdom", "Phylum", "prevalence_frac", "degree", metrics)
stopifnot(all(required %in% names(node.metrics)), all(c("Species", "index") %in% names(euk.indicators)),
          all(c("genome", "index") %in% names(prok.indicators)))

n.original <- nrow(node.metrics)
node.metrics <- node.metrics %>%
  filter(if_all(all_of(required), ~ !is.na(.x))) %>%
  mutate(node = as.character(node), subnetwork = as.character(subnetwork),
         Phylum = str_remove(Phylum, "^p__"), Kingdom = str_remove(Kingdom, "^d__")) %>%
  arrange(subnetwork, node)
message("Retained ", nrow(node.metrics), " / ", n.original, " node records with complete analysis fields")
if (nrow(node.metrics) == 0) stop("No complete node records remain.")
if (anyDuplicated(node.metrics[c("subnetwork", "node")])) stop("Duplicate node/subnetwork records in input.")

conditions <- c("Early", "Late", "Dry", "Wet", "Poor", "Rich")
condition.bits <- set_names(as.integer(2^(seq_along(conditions) - 1)), conditions)
indicators <- bind_rows(euk.indicators %>% transmute(node = str_replace_all(Species, fixed(" "), "."), index),
                        prok.indicators %>% transmute(node = genome, index)) %>%
  filter(index %in% conditions) %>% distinct(node, index)
indicator.profiles <- indicators %>%
  group_by(node) %>% summarise(profile = sum(condition.bits[index]), .groups = "drop")

network.conditions <- tibble(subnetwork = sort(unique(node.metrics$subnetwork))) %>%
  mutate(succession = str_extract(subnetwork, "^(Early|Late)"),
         moisture = str_extract(subnetwork, "(Dry|Wet)"),
         fertility = str_extract(subnetwork, "(Poor|Rich)$"))
if (anyNA(network.conditions[c("succession", "moisture", "fertility")])) stop("Unexpected subnetwork naming convention.")
network.conditions <- network.conditions %>%
  mutate(network.mask = condition.bits[succession] + condition.bits[moisture] + condition.bits[fertility]) %>%
  select(subnetwork, network.mask)

node.metrics <- node.metrics %>%
  left_join(indicator.profiles, by = "node") %>% left_join(network.conditions, by = "subnetwork") %>%
  mutate(profile = replace_na(as.integer(profile), 0L),
         indicator.match = bitwAnd(profile, network.mask) > 0L,
         indicator.any = profile > 0L,
         indicator.class = case_when(indicator.match ~ "Matching", indicator.any ~ "Nonmatching", TRUE ~ "Never"))

taxon.lookup <- node.metrics %>%
  group_by(node) %>%
  summarise(Kingdom = first(Kingdom), Phylum = first(Phylum), n.networks = n_distinct(subnetwork),
            mean.prevalence = mean(prevalence_frac), profile = first(profile), .groups = "drop") %>%
  arrange(node) %>% mutate(taxon.id = row_number())
if (any(node.metrics %>% group_by(node) %>% summarise(n.profile = n_distinct(profile), .groups = "drop") %>% pull(n.profile) != 1))
  stop("Inconsistent indicator profiles for repeated taxa.")
node.metrics <- node.metrics %>% left_join(taxon.lookup %>% select(node, taxon.id), by = "node")

### Within-network metric ranks and nuisance adjustment
network.levels <- sort(unique(node.metrics$subnetwork))
node.metrics <- node.metrics %>%
  mutate(subnetwork = factor(subnetwork, levels = network.levels)) %>%
  group_by(subnetwork) %>%
  mutate(across(all_of(metrics), ~ rank(.x, ties.method = "average") / n(), .names = "rank.{.col}")) %>%
  ungroup()

# Adjust ranks for network-specific quadratic prevalence trends.
# Only betweenness is additionally adjusted for log(1 + total degree).
network.model <- function(dat, metric) {
  rhs <- "0 + subnetwork + subnetwork:prevalence_frac + subnetwork:I(prevalence_frac^2)"
  if (metric == "betweenness") rhs <- paste(rhs, "+ subnetwork:log1p(degree) + subnetwork:I(log1p(degree)^2)")
  residuals(lm(as.formula(paste0("rank.", metric, " ~ ", rhs)), data = dat))
}
for (metric in metrics) node.metrics[[paste0("resid.", metric)]] <- network.model(node.metrics, metric)

### Efficient permutation inputs
# Prepared numeric arrays avoid thousands of repeated dplyr group/pivot operations.
profile <- taxon.lookup$profile
node.taxon <- node.metrics$taxon.id
node.mask <- as.integer(node.metrics$network.mask)
node.network <- as.integer(node.metrics$subnetwork)
network.rows <- split(seq_len(nrow(node.metrics)), node.network)
residuals.mat <- as.matrix(node.metrics[paste0("resid.", metrics)])
min.group <- 5L

score.network <- function(q, contrast) {
  matching <- bitwAnd(q[node.taxon], node.mask) > 0L
  other <- if (contrast == "never") q[node.taxon] == 0L else q[node.taxon] > 0L & !matching
  effects <- matrix(NA_real_, nrow = length(network.rows), ncol = length(metrics))
  for (i in seq_along(network.rows)) {
    ix <- network.rows[[i]]
    a <- ix[matching[ix]]; b <- ix[other[ix]]
    if (length(a) >= min.group && length(b) >= min.group)
      effects[i, ] <- colMeans(residuals.mat[a, , drop = FALSE]) - colMeans(residuals.mat[b, , drop = FALSE])
  }
  dimnames(effects) <- list(network.levels, metrics)
  effects
}

score.profile <- function(q, contrast) {
  effects <- score.network(q, contrast)
  out <- colMeans(effects, na.rm = TRUE)
  out[colSums(!is.na(effects)) == 0] <- NA_real_
  out
}

# The taxon (not its separate network occurrences) is the permutation unit.
# Permutation strata = Phylum x number of observed networks x prevalence bin.
make.strata <- function(taxa, n.bins) {
  taxa %>% group_by(Phylum, n.networks) %>%
    mutate(bin = if (n() >= 3 * n.bins) ntile(rank(mean.prevalence, ties.method = "first"), n.bins) else 1L) %>%
    ungroup() %>% mutate(stratum = paste(Phylum, n.networks, bin, sep = "|"))
}

make.groups <- function(taxa, contrast) {
  taxa %>% filter(contrast == "never" | profile > 0L) %>%
    group_by(stratum) %>% summarise(ids = list(taxon.id), .groups = "drop") %>%
    filter(lengths(ids) > 1L) %>% pull(ids)
}

permute.profiles <- function(q, groups) {
  out <- q
  for (ids in groups) out[ids] <- sample(q[ids], length(ids), replace = FALSE)
  out
}

run.chunk <- function(n, q, groups, contrast) {
  pooled <- matrix(NA_real_, nrow = n, ncol = length(metrics))
  network <- array(NA_real_, dim = c(n, length(network.levels), length(metrics)))
  for (i in seq_len(n)) {
    effects <- score.network(permute.profiles(q, groups), contrast)
    network[i, , ] <- effects
    pooled[i, ] <- colMeans(effects, na.rm = TRUE)
    pooled[i, colSums(!is.na(effects)) == 0] <- NA_real_
  }
  list(pooled = pooled, network = network)
}

run.null <- function(n.bins, contrast, n.perm) {
  taxa <- make.strata(taxon.lookup, n.bins)
  groups <- make.groups(taxa, contrast)
  observed.network <- score.network(profile, contrast)
  observed <- score.profile(profile, contrast)
  if (anyNA(observed)) stop("No usable subnetworks for one or more metrics in: ", contrast)
  sizes <- rep(chunk.size, n.perm %/% chunk.size)
  if (n.perm %% chunk.size) sizes <- c(sizes, n.perm %% chunk.size)
  message("Running: ", n.bins, " bins / ", contrast, " / ", n.perm, " permutations / ", length(groups), " groups")

  # future_lapply supplies reproducible, independent random streams to every chunk.
  # Workers are separate R processes (multisession), including on Windows.
  if (n.cores > 1L) {
    null.list <- future.apply::future_lapply(sizes, run.chunk, q = profile, groups = groups,
                                            contrast = contrast, future.seed = seed + n.bins + (contrast == "nonmatching"))
  } else {
    set.seed(seed + n.bins + (contrast == "nonmatching"))
    null.list <- lapply(sizes, run.chunk, q = profile, groups = groups, contrast = contrast)
  }
  null <- do.call(rbind, lapply(null.list, `[[`, "pooled"))
  network.null <- array(NA_real_, dim = c(n.perm, length(network.levels), length(metrics)))
  offset <- 0L
  for (part in null.list) {
    ix <- seq.int(offset + 1L, offset + dim(part$network)[1])
    network.null[ix, , ] <- part$network
    offset <- max(ix)
  }
  center <- colMeans(null)
  difference <- observed - center
  p <- (1 + colSums(abs(sweep(null, 2, center, "-")) >= rep(abs(difference), each = n.perm))) / (n.perm + 1)
  specification <- paste0("phylum_prevalence_", if_else(n.bins == 3L, "tertiles", "quintiles"))
  results <- tibble(specification, contrast, metric = metrics, observed_effect = observed,
                    null_mean = center, excess_vs_null = difference,
                    null_low = apply(null, 2, quantile, probs = 0.025),
                    null_high = apply(null, 2, quantile, probs = 0.975),
                    p, p_bh = p.adjust(p, method = "BH"), permutations = n.perm)
  coverage <- tibble(specification, contrast, n_taxa = nrow(taxon.lookup),
                     n_indicator_taxa = sum(profile > 0), n_swappable_taxa = length(unique(unlist(groups))),
                     n_swappable_indicator_taxa = sum(profile[unique(unlist(groups))] > 0),
                     n_exchange_groups = length(groups))
  draws <- as_tibble(null, .name_repair = "minimal") %>%
    set_names(metrics) %>% mutate(permutation = row_number()) %>%
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "null_effect") %>%
    mutate(specification = specification, contrast = contrast) %>%
    select(specification, contrast, permutation, metric, null_effect)
  network.draws <- as_tibble(as.data.frame.table(network.null, responseName = "null_effect")) %>%
    transmute(permutation = as.integer(Var1), subnetwork = network.levels[as.integer(Var2)],
              metric = metrics[as.integer(Var3)], null_effect,
              specification = specification, contrast = contrast) %>%
    left_join(as_tibble(as.data.frame.table(observed.network, responseName = "observed_effect")) %>%
                transmute(subnetwork = as.character(Var1), metric = as.character(Var2), observed_effect),
              by = c("subnetwork", "metric")) %>%
    select(specification, contrast, permutation, subnetwork, metric, null_effect, observed_effect)
  list(results = results, coverage = coverage, draws = draws, network.draws = network.draws)
}

### Null models: H1 indicator identity, H2 compatibility, H3 degree-adjusted betweenness
if (n.cores > 1L) {
  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE))
    stop("Install 'future' and 'future.apply' for parallel processing, or set n.cores <- 1.")
  future::plan(future::multisession, workers = n.cores)
}

runs <- list()
for (bins in specifications) for (contrast in c("never", "nonmatching")) {
  key <- paste(bins, contrast, sep = "_")
  runs[[key]] <- run.null(bins, contrast, n.perm)
  readr::write_csv(runs[[key]]$results, file.path(out.dir, paste0("null_", key, ".csv")))
}
if (n.cores > 1L) future::plan(future::sequential)

final.results <- bind_rows(map(runs, "results"))
coverage <- bind_rows(map(runs, "coverage"))
null.draws <- bind_rows(map(runs, "draws"))
network.null.draws <- bind_rows(map(runs, "network.draws"))
readr::write_csv(null.draws, file.path(out.dir, "indicator_final_null_draws.csv"))
readr::write_csv(network.null.draws, file.path(out.dir, "indicator_network_null_draws.csv"))
readr::write_csv(final.results, file.path(out.dir, "indicator_final_null_model_results.csv"))
readr::write_csv(coverage, file.path(out.dir, "indicator_final_null_model_coverage.csv"))

### Descriptive contrasts by kingdom and subnetwork
kingdom.effects <- node.metrics %>%
  filter(indicator.class %in% c("Matching", "Nonmatching")) %>%
  group_by(Kingdom, subnetwork, indicator.class) %>%
  summarise(n = n(), across(starts_with("resid."), mean), .groups = "drop") %>%
  pivot_longer(starts_with("resid."), names_to = "metric", values_to = "mean.residual") %>%
  pivot_wider(names_from = indicator.class, values_from = c(n, mean.residual)) %>%
  mutate(effect = if_else(n_Matching >= min.group & n_Nonmatching >= min.group,
                          mean.residual_Matching - mean.residual_Nonmatching, NA_real_))
readr::write_csv(kingdom.effects, file.path(out.dir, "indicator_kingdom_network_diagnostics.csv"))

indicator.counts <- node.metrics %>% count(subnetwork, Kingdom, indicator.class, name = "n")
readr::write_csv(indicator.counts, file.path(out.dir, "indicator_class_counts.csv"))
message("Done. Results saved to: ", out.dir)
