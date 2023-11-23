# Description of the script

library(optparse)
library(tidyverse)
library(openxlsx)

combined.dd <- read.xlsx("/home/fdelogu/biowide_net/results/split_nets/combined_nets.xlsx", colNames = T)
combined.tax <- rbind((read.csv("/home/fdelogu/biowide_net/data/2023-09-25_plant_tax.csv", header = T) %>%
                         mutate(OTU = X) %>%
                         select(-X)),
                      (read.csv("/home/fdelogu/biowide_net/data/2023-09-25_genus_tax.csv", header = T) %>%
                         mutate(OTU = Genus,
                                Species = NA) %>%
                         select(-X) %>%
                         distinct()))

OTU.split.bin <- rbind(data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/jaccard/EarlyDryPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyDryPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/jaccard/EarlyDryRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyDryRich"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/jaccard/EarlyWetPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyWetPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/jaccard/EarlyWetRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyWetRich"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/jaccard/LateDryPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateDryPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/jaccard/LateDryRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateDryRich"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/jaccard/LateWetPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateWetPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/jaccard/LateWetRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateWetRich"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/sparcc/EarlyDryPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyDryPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/sparcc/EarlyDryRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyDryRich"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/sparcc/EarlyWetPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyWetPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/sparcc/EarlyWetRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyWetRich"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/sparcc/LateDryPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateDryPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/sparcc/LateDryRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateDryRich"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/sparcc/LateWetPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateWetPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/split_OTUs/sparcc/LateWetRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateWetRich"))

OTU.split.tax <- left_join(OTU.split.bin,
                           combined.tax,
                           by="OTU") %>%
  group_by(Kingdom, label) %>%
  summarise(n=n()) %>%
  filter(!is.na(Kingdom)) %>% 
  pivot_wider(names_from = label, id_cols = Kingdom, values_from = n) %>%
  column_to_rownames("Kingdom") %>%
  as.matrix()
OTU.split.tax

combined.dd.tax <- combined.dd %>%
  left_join((combined.tax %>%
               mutate(Kingdom.1 = Kingdom) %>%
               select(OTU, Kingdom.1)),
            by=c("OTU.1" = "OTU")) %>%
  filter(!is.na(Kingdom.1)) %>%
  left_join((combined.tax %>%
               mutate(Kingdom.2 = Kingdom) %>%
               select(OTU, Kingdom.2)),
            by=c("OTU.2" = "OTU")) %>%
  filter(!is.na(Kingdom.2))

combined.dd %>%
  ggplot(aes(x=p.value, y=stat, color=tech)) +
  geom_point(alpha=0.1)

combined.dd %>%
  #filter(abs(stat)>0.5) %>%
  group_by(tech, label) %>%
  summarise(n=n()) %>%
  ungroup() %>%
  pivot_wider(names_from = tech, id_cols = label, values_from = n)

kingdom.dd <- combined.tax %>%
  group_by(Kingdom) %>%
  summarise(n=n())

kingdom.counts <- kingdom.dd$n
names(kingdom.counts) <- kingdom.dd$Kingdom

dd <- combined.dd.tax %>%
  group_by(Kingdom.1, Kingdom.2, tech, label) %>%
  summarise(edges=n()) %>%
  rowwise() %>%
  mutate(nodes = if_else(Kingdom.1 == Kingdom.2,
                         OTU.split.tax[Kingdom.1, label],
                         (OTU.split.tax[Kingdom.1, label]+OTU.split.tax[Kingdom.2, label])),
         potential.edges = if_else(Kingdom.1 == Kingdom.2,
                                   ((OTU.split.tax[Kingdom.1, label]**2-OTU.split.tax[Kingdom.1, label])/2),
                                   OTU.split.tax[Kingdom.1, label]*OTU.split.tax[Kingdom.2, label]),
         frac.edges = edges/potential.edges)

dd %>%
  mutate(net = paste0(Kingdom.1, ":", Kingdom.2)) %>%
  ggplot(aes(x = net, y = nodes)) +
  geom_bar(stat = "identity") +
  facet_grid(label~.)
