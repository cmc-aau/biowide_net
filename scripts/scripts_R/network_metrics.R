#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(igraph)
})

option_list <- list(
  make_option(c("-i", "--input"), dest = "input", type = "character", help = "Input Excel workbook (.xlsx/.xls)."),
  make_option(c("--sheet"), dest = "sheet", type = "character", default = "1",
              help = "Worksheet name or 1-based worksheet index [default: %default]."),
  make_option(c("--source-col"), dest = "source_col", type = "character", default = "source",
              help = "Source-node column [default: %default]."),
  make_option(c("--target-col"), dest = "target_col", type = "character", default = "target",
              help = "Target-node column [default: %default]."),
  make_option(c("--weight-col"), dest = "weight_col", type = "character", default = "weight",
              help = "Signed association-weight column, expected in [-1, 1] [default: %default]."),
  make_option(c("--subnetwork-col"), dest = "subnetwork_col", type = "character", default = "subnetwork",
              help = "Column defining the subnetworks [default: %default]."),
  make_option(c("--node-output"), dest = "node_output", type = "character", default = "node_metrics.csv",
              help = "Node-metric output (.csv, .tsv, or .xlsx) [default: %default]."),
  make_option(c("--edge-output"), dest = "edge_output", type = "character", default = "edge_metrics.csv",
              help = "Edge-metric output (.csv, .tsv, or .xlsx) [default: %default]."),
  make_option(c("--network-output"), dest = "network_output", type = "character", default = "network_metrics.csv",
              help = "Subnetwork-metric output (.csv, .tsv, or .xlsx) [default: %default].")
)

opt <- parse_args(OptionParser(option_list = option_list))

stop_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

write_table <- function(x, path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    write.csv(x, path, row.names = FALSE, na = "")
  } else if (ext %in% c("tsv", "txt")) {
    write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  } else if (ext == "xlsx") {
    if (!requireNamespace("writexl", quietly = TRUE)) {
      stop("Writing .xlsx files requires package 'writexl'. Install it with install.packages('writexl').",
           call. = FALSE)
    }
    writexl::write_xlsx(x, path)
  } else {
    stop("Unsupported output extension for '", path, "'. Use .csv, .tsv, .txt, or .xlsx.",
         call. = FALSE)
  }
}

as_sheet <- function(x) {
  if (grepl("^[0-9]+$", x)) as.integer(x) else x
}

safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

make_graph <- function(edges, source_col, target_col, weight_col) {
  graph_df <- data.frame(
    source = as.character(edges[[source_col]]),
    target = as.character(edges[[target_col]]),
    weight = as.numeric(edges[[weight_col]]),
    edge_row_id = seq_len(nrow(edges)),
    stringsAsFactors = FALSE
  )
  graph_from_data_frame(graph_df, directed = FALSE)
}

get_node_metrics <- function(g, subnetwork_id) {
  pos_e <- E(g)[weight > 0]
  neg_e <- E(g)[weight < 0]

  g_pos <- subgraph.edges(g, pos_e, delete.vertices = FALSE)
  g_neg <- subgraph.edges(g, neg_e, delete.vertices = FALSE)

  degree_total <- degree(g)
  degree_pos <- degree(g_pos)
  degree_neg <- degree(g_neg)

  strength_pos <- strength(g_pos, weights = E(g_pos)$weight)
  strength_neg <- if (ecount(g_neg) > 0) {
    strength(g_neg, weights = abs(E(g_neg)$weight))
  } else {
    setNames(rep(0, vcount(g_neg)), V(g_neg)$name)
  }

  btw <- betweenness(g, directed = FALSE, weights = NA, normalized = TRUE)
  clustering_pos <- transitivity(g_pos, type = "local", isolates = "NaN")
  clustering_pos[is.nan(clustering_pos)] <- NA_real_

  data.frame(
    subnetwork = subnetwork_id,
    node = V(g)$name,
    degree = as.numeric(degree_total[V(g)$name]),
    positive_degree = as.numeric(degree_pos[V(g)$name]),
    negative_degree = as.numeric(degree_neg[V(g)$name]),
    positive_strength = as.numeric(strength_pos[V(g)$name]),
    negative_strength = as.numeric(strength_neg[V(g)$name]),
    betweenness = as.numeric(btw[V(g)$name]),
    positive_clustering = as.numeric(clustering_pos[V(g)$name]),
    stringsAsFactors = FALSE
  )
}

get_edge_metrics <- function(g, subnetwork_id, original_edges) {
  eb <- edge_betweenness(g, directed = FALSE, weights = NA)
  n <- vcount(g)
  denom <- n * (n - 1) / 2
  eb_norm <- if (denom > 0) eb / denom else rep(NA_real_, length(eb))

  edge_ends <- ends(g, E(g), names = TRUE)
  edge_rows <- as.integer(E(g)$edge_row_id)

  out <- original_edges[edge_rows, , drop = FALSE]
  out$subnetwork <- subnetwork_id
  out$source_node <- edge_ends[, 1]
  out$target_node <- edge_ends[, 2]
  out$weight <- E(g)$weight
  out$sign <- ifelse(E(g)$weight > 0, "positive",
                     ifelse(E(g)$weight < 0, "negative", "zero"))
  out$edge_betweenness <- as.numeric(eb)
  out$edge_betweenness_norm <- as.numeric(eb_norm)

  preferred <- c("subnetwork", "source_node", "target_node", "weight", "sign",
                 "edge_betweenness", "edge_betweenness_norm")
  out[, c(preferred, setdiff(names(out), preferred)), drop = FALSE]
}

get_network_metrics <- function(g, subnetwork_id) {
  pos_e <- E(g)[weight > 0]
  neg_e <- E(g)[weight < 0]
  g_pos <- subgraph.edges(g, pos_e, delete.vertices = FALSE)
  g_neg <- subgraph.edges(g, neg_e, delete.vertices = FALSE)

  n_nodes <- vcount(g)
  n_edges <- ecount(g)
  n_pos <- length(pos_e)
  n_neg <- length(neg_e)

  deg <- degree(g)
  pos_strength <- strength(g_pos, weights = E(g_pos)$weight)
  neg_strength <- if (ecount(g_neg) > 0) {
    strength(g_neg, weights = abs(E(g_neg)$weight))
  } else {
    rep(0, vcount(g_neg))
  }

  trans_pos <- transitivity(g_pos, type = "global")
  if (is.nan(trans_pos)) trans_pos <- NA_real_

  positive_modularity <- NA_real_
  n_positive_communities <- NA_integer_
  if (ecount(g_pos) > 0) {
    non_isolated <- V(g_pos)[degree(g_pos) > 0]
    if (length(non_isolated) >= 2) {
      g_pos_active <- induced_subgraph(g_pos, non_isolated)
      comm <- cluster_louvain(g_pos_active, weights = E(g_pos_active)$weight)
      positive_modularity <- modularity(comm)
      n_positive_communities <- length(unique(membership(comm)))
    }
  }

  comps <- components(g)

  data.frame(
    subnetwork = subnetwork_id,
    n_nodes = n_nodes,
    n_edges = n_edges,
    n_positive_edges = n_pos,
    n_negative_edges = n_neg,
    n_zero_edges = n_edges - n_pos - n_neg,
    positive_edge_fraction = if (n_edges > 0) n_pos / n_edges else NA_real_,
    density = edge_density(g, loops = FALSE),
    mean_degree = safe_mean(deg),
    mean_positive_strength = safe_mean(pos_strength),
    mean_negative_strength = safe_mean(neg_strength),
    positive_transitivity = trans_pos,
    n_components = comps$no,
    largest_component_fraction = if (n_nodes > 0) max(comps$csize) / n_nodes else NA_real_,
    positive_modularity = positive_modularity,
    n_positive_communities = n_positive_communities,
    stringsAsFactors = FALSE
  )
}

analyze_subnetwork <- function(edges, subnetwork_id, source_col, target_col, weight_col) {
  g <- make_graph(edges, source_col, target_col, weight_col)
  list(
    nodes = get_node_metrics(g, subnetwork_id),
    edges = get_edge_metrics(g, subnetwork_id, edges),
    network = get_network_metrics(g, subnetwork_id)
  )
}

stop_if(is.null(opt$input), "Missing required argument: --input")
stop_if(!file.exists(opt$input), paste0("Input file does not exist: ", opt$input))

message("Reading: ", opt$input)
dat <- read_excel(opt$input, sheet = as_sheet(opt$sheet))
dat <- as.data.frame(dat, stringsAsFactors = FALSE)

message("Using columns: source='", opt$source_col, "', target='", opt$target_col,
        "', weight='", opt$weight_col, "', subnetwork='", opt$subnetwork_col, "'")

required_cols <- c(opt$source_col, opt$target_col, opt$weight_col, opt$subnetwork_col)
missing_cols <- setdiff(required_cols, names(dat))
stop_if(length(missing_cols) > 0,
        paste0("Missing required column(s): ", paste(missing_cols, collapse = ", ")))

dat[[opt$source_col]] <- trimws(as.character(dat[[opt$source_col]]))
dat[[opt$target_col]] <- trimws(as.character(dat[[opt$target_col]]))
dat[[opt$subnetwork_col]] <- as.character(dat[[opt$subnetwork_col]])
dat[[opt$weight_col]] <- suppressWarnings(as.numeric(dat[[opt$weight_col]]))

bad_missing <- is.na(dat[[opt$source_col]]) | dat[[opt$source_col]] == "" |
               is.na(dat[[opt$target_col]]) | dat[[opt$target_col]] == "" |
               is.na(dat[[opt$subnetwork_col]]) | dat[[opt$subnetwork_col]] == "" |
               is.na(dat[[opt$weight_col]])
stop_if(any(bad_missing),
        paste0("Found ", sum(bad_missing),
               " row(s) with missing/invalid source, target, subnetwork, or weight values."))

stop_if(any(dat[[opt$weight_col]] < -1 | dat[[opt$weight_col]] > 1),
        "Weights outside [-1, 1] were found.")

self_loops <- dat[[opt$source_col]] == dat[[opt$target_col]]
stop_if(any(self_loops),
        paste0("Found ", sum(self_loops),
               " self-loop(s). Remove them before running this script."))

dup_key <- paste(
  dat[[opt$subnetwork_col]],
  pmin(dat[[opt$source_col]], dat[[opt$target_col]]),
  pmax(dat[[opt$source_col]], dat[[opt$target_col]]),
  sep = "\r"
)
stop_if(anyDuplicated(dup_key) > 0,
        "Duplicate undirected edges were found within at least one subnetwork. Resolve duplicates before analysis.")

if (any(dat[[opt$weight_col]] == 0)) {
  warning("Zero-weight edges are retained in the full topology but excluded from positive/negative strength and positive-network metrics.",
          call. = FALSE)
}

subnetwork_ids <- unique(dat[[opt$subnetwork_col]])
message("Found ", length(subnetwork_ids), " subnetwork(s): ",
        paste(subnetwork_ids, collapse = ", "))

results <- lapply(subnetwork_ids, function(id) {
  edges <- dat[dat[[opt$subnetwork_col]] == id, , drop = FALSE]
  message("Analyzing subnetwork '", id, "' (", nrow(edges), " edges)")
  analyze_subnetwork(edges, id, opt$source_col, opt$target_col, opt$weight_col)
})

node_metrics <- do.call(rbind, lapply(results, `[[`, "nodes"))
edge_metrics <- do.call(rbind, lapply(results, `[[`, "edges"))
network_metrics <- do.call(rbind, lapply(results, `[[`, "network"))

row.names(node_metrics) <- NULL
row.names(edge_metrics) <- NULL
row.names(network_metrics) <- NULL

write_table(node_metrics, opt$node_output)
write_table(edge_metrics, opt$edge_output)
write_table(network_metrics, opt$network_output)

message("Done.")
message("Node metrics:    ", opt$node_output)
message("Edge metrics:    ", opt$edge_output)
message("Network metrics: ", opt$network_output)
