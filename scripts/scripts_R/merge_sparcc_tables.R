# Pearson's correlation for continuous and dichotomous variables, equivalent to point biserial correlation
# Used to link the an OTU table with counts and one with presence/absence values

library(optparse)
library(tidyverse)
library(openxlsx)

option_list = list(
  make_option(c("-r", "--cor_mat"), type="character", default=NULL, 
              help="Correlation table from sparcc (.csv).", metavar="character"),
  make_option(c("-v", "--cov_mat"), type="character", default=NULL, 
              help="Covariance table from sparcc (.csv).", metavar="character"),
  make_option(c("-s", "--one_sided_p_mat"), type="character", default=NULL, 
              help="One sided p value table from sparcc (.csv).", metavar="character"),
  make_option(c("-t", "--two_sided_p_mat"), type="character", default=NULL, 
              help="Two sided p value table from sparcc (.csv).", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", 
              help="Output for the combined table.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if(file.exists(opt$out_file)){
  file.remove(opt$out_file)
}

wb <- createWorkbook()

addWorksheet(wb, "correlation")
writeData(wb, "correlation", read.csv(opt$cor_mat, row.names = 1), rowNames = T)

#addWorksheet(wb, "covariance")
#writeData(wb, "covariance", read.csv(opt$cov_mat, row.names = 1), rowNames = T)

addWorksheet(wb, "1s.p.value")
writeData(wb, "1s.p.value", read.csv(opt$one_sided_p_mat, row.names = 1), rowNames = T)

#addWorksheet(wb, "2s.p.value")
#writeData(wb, "2s.p.value", read.csv(opt$two_sided_p_mat, row.names = 1), rowNames = T)

saveWorkbook(wb, opt$out_file)
