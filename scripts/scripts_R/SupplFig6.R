# Panel B: eight subnetworks as boxplots with observed values and significance stars

### Loading libraries
suppressMessages(library(tidyverse))

### Settings
base.dir <- "/home/fdelogu/biowide_net"
out.dir <- file.path(base.dir, "results", "indicator_null_models")
fig.dir <- file.path(base.dir, "results", "figures")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

specification <- "phylum_prevalence_tertiles"
adjust.p <- TRUE
null.file <- file.path(out.dir, "indicator_network_null_draws.csv")

network.levels <- c("EarlyDryPoor", "EarlyDryRich", "EarlyWetPoor", "EarlyWetRich",
                    "LateDryPoor", "LateDryRich", "LateWetPoor", "LateWetRich")
network.labels <- c("EarlyDryPoor", "EarlyDryRich", "EarlyWetPoor", "EarlyWetRich",
                    "LateDryPoor", "LateDryRich", "LateWetPoor", "LateWetRich")
network.colors <- set_names(c("#E7D9A5", "#B8D7B0", "#B8D9E8", "#F6D0A5",
                              "#98803B", "#34734C", "#246A8A", "#C4772B"), network.levels)
metric.labels <- c("positive_degree" = "Positive degree",
                   "positive_strength" = "Positive strength",
                   "betweenness" = "Degree-adjusted betweenness")
contrast.labels <- c("never" = "Matching vs. never indicators",
                     "nonmatching" = "Matching vs. nonmatching indicators")

### Loading network-specific null distributions
panel.b.data <- read_csv(null.file, show_col_types = FALSE) %>%
  filter(.data$specification == .env$specification,
         contrast %in% names(contrast.labels), metric %in% names(metric.labels)) %>%
  mutate(subnetwork = factor(subnetwork, levels = network.levels, labels = network.labels),
         contrast = factor(contrast, levels = names(contrast.labels), labels = contrast.labels),
         metric = factor(metric, levels = names(metric.labels), labels = metric.labels))

stopifnot(nrow(panel.b.data) > 0,
          all(c("permutation", "null_effect", "observed_effect") %in% names(panel.b.data)))

observed <- panel.b.data %>%
  distinct(contrast, metric, subnetwork, observed_effect) %>%
  filter(!is.na(observed_effect))

ecospace.order <- c("EarlyDryPoor", "LateDryPoor",
                    "EarlyWetPoor", "LateWetPoor",
                    "EarlyDryRich", "LateDryRich",
                    "EarlyWetRich", "LateWetRich")

### Subnetwork-specific permutation p-values and significance stars
sig.table <- panel.b.data %>%
  filter(!is.na(null_effect), !is.na(observed_effect)) %>%
  group_by(contrast, metric, subnetwork) %>%
  summarise(observed_effect = first(observed_effect),
            null.center = mean(null_effect),
            n.perm = n(),
            p.value = (sum(abs(null_effect - mean(null_effect)) >= abs(first(observed_effect) - mean(null_effect))) + 1) / (n() + 1),
            q05 = quantile(null_effect, 0.05),
            q95 = quantile(null_effect, 0.95),
            null.min = min(null_effect),
            null.max = max(null_effect),
            .groups = "drop") %>%
  group_by(contrast, metric) %>%
  mutate(p.plot = if (adjust.p) p.adjust(p.value, method = "BH") else p.value,
         star = case_when(p.plot <= 0.001 ~ "***",
                          p.plot <= 0.01 ~ "**",
                          p.plot <= 0.05 ~ "*",
                          TRUE ~ "")) %>%
  ungroup() %>%
  mutate(range.all = pmax(abs(null.min), abs(null.max), abs(observed_effect)),
         star.y = if_else(observed_effect >= null.center,
                          pmax(null.max, observed_effect) + 0.08 * range.all,
                          pmin(null.min, observed_effect) - 0.10 * range.all),
         vjust = if_else(observed_effect >= null.center, 0, 1))

panel.b <- panel.b.data %>%
  mutate(subnetwork = factor(subnetwork, levels = ecospace.order)) %>%
  filter(!is.na(null_effect)) %>%
  ggplot(aes(x = subnetwork, y = null_effect, fill = subnetwork, color = subnetwork)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dotted", color = "grey55") +
  geom_boxplot(width = 0.7, linewidth = 0.35, alpha = 0.45, outlier.shape = NA) +
  geom_point(data = observed,
             aes(x = subnetwork, y = observed_effect, color = subnetwork),
             inherit.aes = FALSE, shape = 4, size = 2.8, stroke = 1.0) +
  geom_text(data = sig.table %>% filter(star != ""),
            aes(x = subnetwork, y = star.y, label = star, color = subnetwork, vjust = vjust),
            inherit.aes = FALSE, size = 4, fontface = "bold", show.legend = FALSE) +
  facet_grid(contrast ~ metric, scales = "free_y") +
  scale_fill_manual(values = network.colors, drop = FALSE) +
  scale_color_manual(values = network.colors, drop = FALSE) +
  labs(x = "Ecospace",
       y = "Prevalence-adjusted mean percentile-rank difference") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

print(panel.b)

ggsave(file.path(fig.dir, "SupplFig7.png"), panel.b,
       dpi = 300, width = 24, height = 17, units = "cm", bg = "white")
ggsave(file.path(fig.dir, "SupplFig7.svg"), panel.b,
       width = 24, height = 17, units = "cm", bg = "white")

write_csv(sig.table, file.path(fig.dir, "SupplFig7_pval.csv"))
