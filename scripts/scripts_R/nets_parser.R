# Script description

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(openxlsx))

option_list = list(
  make_option(c("-i", "--input_files"), type="character", 
              help="Comma separated list of net files.", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", 
              help="Output for the combined table.", metavar="character"),
  make_option(c("-p", "--p_threshold"), type="double",
              default = 10,
              help="P value threshold for filtering.", metavar="double"),
  make_option(c("-s", "--stat_threshold"), type="double", 
              default = -1000,
              help="Stat value threshold for filtering.", metavar="double")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

net.table.list <- unlist(str_split(opt$input_files, ","))
#print(net.table.list)

#net.table.list <- c("/home/fdelogu/biowide_net/results/split_nets/jaccard/EarlyDryPoor_final.xlsx", "/home/fdelogu/biowide_net/results/split_nets/jaccard/EarlyWetPoor_final.xlsx")

net.reformat <- function(fname){

  # Reformats the network to a long shape and reduces the size using the symmetry of the net and eliminating the diagonal

  dd.stats <- read.xlsx(fname, colNames = T, rowNames = T, sheet = 1)
  dd.pvalues <- read.xlsx(fname, colNames = T, rowNames = T, sheet = 2)
  
  long.stats <- dd.stats %>%
    rownames_to_column("OTU.A") %>%
    pivot_longer(names_to = "OTU.B", values_to = "stat", -OTU.A)
  
  long.pvalues <- dd.pvalues %>%
    rownames_to_column("OTU.A") %>%
    pivot_longer(names_to = "OTU.B", values_to = "p.value", -OTU.A)
  
  long.combined <- inner_join(long.stats,
                              long.pvalues,
                              by=c("OTU.A", "OTU.B")) %>%
    filter(OTU.A!=OTU.B) %>%
    rowwise() %>%
    mutate(OTU.1 = sort(c(OTU.A, OTU.B))[1],
           OTU.2 = sort(c(OTU.A, OTU.B))[2]) %>%
    ungroup() %>%
    select(-OTU.A, -OTU.B) %>%
    distinct()
  
  return(long.combined)
}

net.filter <- function(dd.long, p.threshold=10, stat.threshold=-10){
  
  dd.filtered <- dd.long %>%
    filter(p.value < p.threshold,
           stat > stat.threshold)
  
  return(dd.filtered)
  
}

net.grouping <- function(dd.long, fname){
  
  tech <- unlist(str_split(fname, "/"))
  tech <- tech[length(tech)-1]
  
  label <- unlist(str_split(fname, "/"))
  label <- label[length(label)]
  label <- unlist(str_split(label, "_"))[1]
  
  dd.return <- dd.long %>%
    mutate(tech = tech,
           label = label)
  
  return(dd.return)
  
}

file2net <- function(fname, p.threshold=10, stat.threshold=-10){
  
  net.long <- net.reformat(fname)
  net.filtered <- net.filter(net.long, p.threshold=p.threshold, stat.threshold=stat.threshold)
  net.annotated <- net.grouping(net.filtered, fname)
  
  return(net.annotated)
  
}

combined.nets <- lapply(net.table.list, file2net, p.threshold=opt$p_threshold, stat.threshold=opt$stat_threshold) %>%
  bind_rows(.id = "column_label")


if(file.exists(opt$out_file)){
  file.remove(opt$out_file)
}

wb <- createWorkbook()

addWorksheet(wb, "combined")
writeData(wb, "combined", combined.nets, rowNames = F)

saveWorkbook(wb, opt$out_file)
