# Split OTU table into subsets according to a mapping file

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(matrixStats))

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="OTU table (.csv)", metavar="character"),
  make_option(c("-m", "--mapping"), type="character", 
              help="Mapping file (.tsv). The first column contains the sample IDs, the second the mapping. The header is required.", metavar="character"),
  make_option(c("-o", "--out_folder"), type="character", 
              help="Output folder.", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

otu.table <- read.csv(opt$file, row.names = 1)

map.table <- read_tsv(opt$mapping) %>%
  `colnames<-`(c("ID", "map_label"))

output_dir <- opt$out_folder

if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = T)
}

#head(otu.table)

#head(map.table)

indices <- unique(map.table %>% pull(map_label))

for(index in indices){

  print(index)
  
  samples.index <- map.table %>%
    filter(map_label==index) %>%
    pull(ID)
  
  otu.index <- otu.table[,samples.index]
  
  print(dim(otu.index))
  
  otu.index.reduced <- otu.index[rowSums(as.matrix(otu.index))>1,]

  print(dim(otu.index.reduced))
 
  write.csv(otu.index.reduced, file=paste0(output_dir, "/", index, "_OTU.csv"))
     
}

