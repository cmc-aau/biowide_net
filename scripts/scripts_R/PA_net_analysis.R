# Analysis of the networks

## Setting

### Loading libraries

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(openxlsx))

##S Loading data

combined.dd <- read.xlsx("/home/fdelogu/biowide_net/results/PA_split_nets/combined_nets.xlsx", colNames = T)
combined.tax <- rbind((read.csv("/home/fdelogu/biowide_net/data/2023-09-25_plant_tax.csv", header = T) %>%
                         mutate(OTU = X) %>%
                         select(-X)),
                      (read.csv("/home/fdelogu/biowide_net/data/2023-09-25_genus_tax.csv", header = T) %>%
                         mutate(OTU = Genus,
                                Species = NA) %>%
                         select(-X) %>%
                         distinct()))

OTU.split.bin <- rbind(data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/PA_split_OTUs/jaccard/EarlyDryPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyDryPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/PA_split_OTUs/jaccard/EarlyDryRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyDryRich"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/PA_split_OTUs/jaccard/EarlyWetPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyWetPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/PA_split_OTUs/jaccard/EarlyWetRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "EarlyWetRich"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/PA_split_OTUs/jaccard/LateDryPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateDryPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/PA_split_OTUs/jaccard/LateDryRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateDryRich"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/PA_split_OTUs/jaccard/LateWetPoor_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateWetPoor"),
                       data.frame(OTU = read.csv("/home/fdelogu/biowide_net/results/PA_split_OTUs/jaccard/LateWetRich_OTU.csv",
                                                 header = T, row.names = 1) %>%
                                    rownames(),
                                  label = "LateWetRich"))


micro.ind <- left_join(read.csv("/home/fdelogu/biowide_net/data/2023-09-25_indicator_genera_rel.csv", header = T),
                       (combined.tax %>%
                          filter(Kingdom %in% c("Archaea", "Bacteria")) %>%
                          select(Genus, Species, OTU)),
                       by="Genus")
above.ind <- left_join(read.csv("/home/fdelogu/biowide_net/data/2023-09-25_indicator_plants.csv", header = T),
          (combined.tax %>%
             filter(!Kingdom %in% c("Archaea", "Bacteria")) %>%
             select(Species, OTU)),
          by="Species")
combined.ind <- rbind(micro.ind, above.ind)


## Exploratory

### Compile list of OTUs per condition per subnetwork

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

#combined.dd %>%
#  ggplot(aes(x=p.value, y=stat, color=tech)) +
#  geom_point(alpha=0.1)

#combined.dd %>%
  #filter(abs(stat)>0.5) %>%
#  group_by(tech, label) %>%
#  summarise(n=n()) %>%
#  ungroup() %>%
#  pivot_wider(names_from = tech, id_cols = label, values_from = n)

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
  ggplot(aes(x = net, y = frac.edges)) +
  geom_bar(stat = "identity") +
  facet_grid(label~.)

### Indicator + network analysis

nodes.dd <- rbind(data.frame(OTU=combined.dd$OTU.1, label=combined.dd$label),
                  data.frame(OTU=combined.dd$OTU.2, label=combined.dd$label))

nodes.ind <- left_join(nodes.dd,
          (combined.ind %>%
             select(OTU, index)),
          by="OTU",
          relationship = "many-to-many") %>%
  filter(!is.na(index)) %>%
  rowwise() %>%
  mutate(ind.node = if_else(str_detect(label, index), 1, 0)) %>%
  distinct()

combined.dd.ind <- combined.dd.tax %>%
  left_join((nodes.ind %>%
               filter(ind.node==1) %>%
               mutate(ind.node.1 = ind.node) %>%
               select(-ind.node, -index)),
            by=c("OTU.1"="OTU",
                 "label"),
            relationship="many-to-many") %>%
  left_join((nodes.ind %>%
               filter(ind.node==1) %>%
               mutate(ind.node.2 = ind.node) %>%
               select(-ind.node, -index)),
            by=c("OTU.2"="OTU",
                 "label"),
            relationship="many-to-many") %>%
  mutate(ind.node.1 = if_else(is.na(ind.node.1), 0, 1),
         ind.node.2 = if_else(is.na(ind.node.2), 0, 1))

combined.dd.ind %>%
  mutate(ind.node.3=ind.node.1+ind.node.2) %>%
  group_by(label, ind.node.3) %>%
  summarise(n=n()) %>%
  pivot_wider(id_cols = label, names_from = ind.node.3, values_from = n)

edges.strata.between <- combined.dd.ind %>%
  mutate(subnet.1 = if_else(Kingdom.1%in%c("Archaea", "Bacteria"), "Below", "Above"),
         subnet.2 = if_else(Kingdom.2%in%c("Archaea", "Bacteria"), "Below", "Above"),
         edge.type = if_else(subnet.1==subnet.2, "intra", "between")) %>%
  mutate(ind.node.3 = ind.node.1 + ind.node.2) %>%
  group_by(label, edge.type, ind.node.3) %>%
  summarise(n=n()) %>%
  ggplot(aes(x=label, y=n, fill=edge.type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(.~as.factor(ind.node.3), scales = "free") +
  coord_flip() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=0.75)) +
  labs(x="Stratum", y="Number of edges")

jac.strata.between <- combined.dd.ind %>%
  mutate(subnet.1 = if_else(Kingdom.1%in%c("Archaea", "Bacteria"), "Below", "Above"),
         subnet.2 = if_else(Kingdom.2%in%c("Archaea", "Bacteria"), "Below", "Above"),
         edge.type = if_else(subnet.1==subnet.2, "intra", "between")) %>%
  mutate(ind.node.3 = ind.node.1 + ind.node.2) %>%
  filter(stat > 0) %>%
  ggplot(aes(x = label, y=stat, fill = edge.type)) +
  geom_violin(position = position_dodge(0.75)) +
  geom_boxplot(position = position_dodge(0.75), width = 0.125, alpha = 0.5) +
  facet_grid(.~as.factor(ind.node.3), scales = "free") +
  coord_flip() +
  theme_classic() +
  labs(x="Stratum", y="Modified Jaccard similarity")

ggsave("/home/fdelogu/biowide_net/results/figures/edges_between.png", edges.strata.between, width = 20, height = 10, units = "cm")
ggsave("/home/fdelogu/biowide_net/results/figures/jac_between.png", jac.strata.between, width = 20, height = 10, units = "cm")

combined.dd.ind %>%
  filter(label=="LateDryRich") %>%
  select(-column_label, -tech) %>%
  write.csv("/home/fdelogu/biowide_net/results/tables/LateDryRich_net.csv")

combined.dd.ind %>%
  filter(label=="LateDryRich",
         stat>0.5 | stat<(-0.5)) %>%
  select(-column_label, -tech) %>%
  write.csv("/home/fdelogu/biowide_net/results/tables/LateDryRich_net_050.csv")

nodes.only <- data.frame(node=c(combined.dd.ind$OTU.1[combined.dd.ind$label=="LateDryRich" & (combined.dd.ind$stat>0.5 | combined.dd.ind$stat<(-0.5))],
                                combined.dd.ind$OTU.2[combined.dd.ind$label=="LateDryRich" & (combined.dd.ind$stat>0.5 | combined.dd.ind$stat<(-0.5))]),
          ind.node=c(combined.dd.ind$ind.node.1[combined.dd.ind$label=="LateDryRich" & (combined.dd.ind$stat>0.5 | combined.dd.ind$stat<(-0.5))],
                     combined.dd.ind$ind.node.2[combined.dd.ind$label=="LateDryRich" & (combined.dd.ind$stat>0.5 | combined.dd.ind$stat<(-0.5))]),
          kingdom=c(combined.dd.ind$Kingdom.1[combined.dd.ind$label=="LateDryRich" & (combined.dd.ind$stat>0.5 | combined.dd.ind$stat<(-0.5))],
                    combined.dd.ind$Kingdom.2[combined.dd.ind$label=="LateDryRich" & (combined.dd.ind$stat>0.5 | combined.dd.ind$stat<(-0.5))])) %>%
  group_by(node, kingdom) %>%
  summarise(indicator.node=if_else(sum(ind.node)>0, 1, 0)) %>%
  ungroup()

nodes.only %>%
  write.csv("/home/fdelogu/biowide_net/results/tables/LateDryRich_nodes_050.csv")
