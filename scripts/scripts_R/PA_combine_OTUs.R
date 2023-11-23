# Merge binary and counts OTU tables into a single binary one

suppressMessages(library(optparse))
suppressMessages(library(tidyverse))
suppressMessages(library(matrixStats))

options(stringsAsFactors = F)


option_list = list(
  make_option(c("-b", "--file_bin"), type="character", default=NULL, 
              help="Binary OTU table (.csv)", metavar="character"),
  make_option(c("-c", "--file_counts"), type="character", default=NULL, 
              help="Counts OTU table (.csv)", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", 
              help="Output combined PA table (.xlsx).", metavar="character"),
  make_option(c("-t", "--pa_threshold"), type="numeric", default=0, 
              help="Threshold for counting presences", metavar="numeric")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

bin.table <- read.csv(opt$file_bin, row.names = 1)
counts.table <- read.csv(opt$file_counts, row.names = 1)

#dim(bin.table)
#dim(counts.table)

combined.table <- rbind(bin.table,
                        counts.table[,colnames(bin.table)]) %>%
  mutate(across(everything(), ~ if_else(.x > opt$pa_threshold, 1, 0)))

#dim(combined.table)
#sum(combined.table!=0&combined.table!=1)

#head(combined.table)

if(file.exists(opt$out_file)){
  file.remove(opt$out_file)
}

write.csv(combined.table, file=opt$out_file)


