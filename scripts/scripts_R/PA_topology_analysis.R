# Analysis of the networks

## Setting

### Loading libraries

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(openxlsx))
suppressMessages(library(ggpubr))

##S Loading data

combined.tax <- rbind((read.csv("/home/fdelogu/biowide_net/data/2023-09-25_plant_tax.csv", header = T) %>%
                         mutate(OTU = X) %>%
                         select(-X)),
                      (read.csv("/home/fdelogu/biowide_net/data/2023-09-25_genus_tax.csv", header = T) %>%
                         mutate(OTU = Genus,
                                Species = NA) %>%
                         select(-X) %>%
                         distinct()))

topology.df <- read.csv("/home/fdelogu/biowide_net/results/tables/LDR_nodes_analysis.csv")

## Exploratory

topology.df %>%
  mutate(name=gsub('"', '', name),
         kingdom=gsub('"', '', kingdom)) %>%
  filter(kingdom!="Archaea") %>%
  select(name,
         kingdom,
         indicator.node,
         AverageShortestPathLength,
         BetweennessCentrality,
         ClosenessCentrality,
         ClusteringCoefficient,
         Degree,
         Eccentricity,
         NeighborhoodConnectivity,
         PartnerOfMultiEdgedNodePairs,
         Radiality,
         Stress,
         TopologicalCoefficient) %>%
  pivot_longer(names_to = "topological.variable",
               values_to = "value",
               -c(name, kingdom, indicator.node)) %>%
  ggplot(aes(x=as.factor(indicator.node), y=value, color=as.factor(indicator.node))) +
  geom_jitter(alpha=0.1) +
  geom_violin(alpha=0.5) +
  geom_boxplot(width = 0.125, alpha = 0.5) +
  stat_compare_means(method = "t.test") +
  facet_wrap(kingdom~topological.variable, scales="free_y", nrow=4) +
  theme_minimal()

## Plot topological metrics

topology.small <- topology.df %>%
  mutate(name=gsub('"', '', name),
         kingdom=gsub('"', '', kingdom),
         `Av. Shortest Path` = AverageShortestPathLength,
         `Clustering Coefficient` = ClusteringCoefficient,
         `Neighborhood Connectivity` = NeighborhoodConnectivity) %>%
  filter(kingdom!="Archaea") %>%
  select(name,
         kingdom,
         indicator.node,
         `Av. Shortest Path`,
         `Clustering Coefficient`,
         Degree,
         `Neighborhood Connectivity`) %>%
  pivot_longer(names_to = "topological.variable",
               values_to = "value",
               -c(name, kingdom, indicator.node)) %>%
  ggplot(aes(x=as.factor(indicator.node), y=value, color=as.factor(indicator.node))) +
  geom_jitter(alpha=0.25) +
  geom_violin(alpha=0.5) +
  geom_boxplot(width = 0.125, alpha = 0.5) +
  stat_compare_means(method = "t.test",
                     label = "p.signif",
                     label.y.npc = 0.4,
                     label.y.npc = 0.95,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.05, Inf), 
                                        symbols = c("p<0.001", "p<0.05", "ns"))) +
  facet_wrap(kingdom~topological.variable, scales="free_y", nrow=4) +
  theme_minimal() +
  theme(legend.position = "bottom")

topology.small
ggsave("/home/fdelogu/biowide_net/results/figures/topology_small.png", topology.small, dpi=300, width = 15, height = 20, units = "cm")
ggsave("/home/fdelogu/biowide_net/results/figures/topology_small.svg", topology.small, width = 20, height = 10, units = "cm")
