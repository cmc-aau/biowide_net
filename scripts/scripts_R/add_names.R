# Add entities' names to sparcc output using the input

library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Sparcc output table (.txt/.csv)", metavar="character"),
  make_option(c("-s", "--names_source"), type="character", 
              help="Source of the names, i.e. the original OTU table provided to sparcc.", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", 
              help="Output file.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

otu.table <- read.csv(opt$names_source, row.names = 1)

to_print <- read.csv(opt$file, row.names = 1) %>%
  `colnames<-`(rownames(otu.table)) %>%
  `rownames<-`(rownames(otu.table))

output_file <- opt$out_file

write.csv(to_print, file=output_file)


