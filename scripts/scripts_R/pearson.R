# Pearson's correlation for continuous and dichotomous variables, equivalent to point biserial correlation
# Used to link the an OTU table with counts and one with presence/absence values

library(optparse)
library(tidyverse)
library(openxlsx)

option_list = list(
  make_option(c("-c", "--counts_OTUs"), type="character", default=NULL, 
              help="OTU table with counts data (.csv).", metavar="character"),
  make_option(c("-b", "--binary_OTUs"), type="character", default=NULL, 
              help="OTU table with presence/absence data (.csv).", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", 
              help="Output file for the correlation similarity with the p value.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

otu.counts <- read.csv(opt$counts_OTUs, row.names = 1)[1:20,]
otu.bin <- read.csv(opt$binary_OTUs, row.names = 1)[1:20,]

## Make the matrices that will store the results (#OTUs in counts table x #OTUs in binary table)

v.mat <- matrix(data = 0,
                nrow = nrow(otu.counts),
                ncol = nrow(otu.bin)) %>%
  `colnames<-`(rownames(otu.bin)) %>%
  `rownames<-`(rownames(otu.counts))

p.mat <- matrix(data = 0,
                nrow = nrow(otu.counts),
                ncol = nrow(otu.bin)) %>%
  `colnames<-`(rownames(otu.bin)) %>%
  `rownames<-`(rownames(otu.counts))


## Compute Pearson correlation and its p value

for(i in 1:(nrow(otu.counts))){
  for(j in 1:nrow(otu.bin)){
    cor.val <- cor.test(as.numeric(otu.counts[i, ]), as.numeric(otu.bin[j, ]), method = "pearson")
    #print(c(rownames(otu.counts)[i], rownames(otu.bin)[j]))
    #print(c(sum(is.na(otu.table[i, ])), sum(is.na(otu.table[j, ]))))
    #print(cor.val)
    v.mat[i, j] <- cor.val$estimate
    p.mat[i, j] <- cor.val$p.value
  }
}


if(file.exists(opt$out_file)){
  file.remove(opt$out_file)
}

wb <- createWorkbook()

addWorksheet(wb, "pearson")
writeData(wb, "pearson", v.mat, rowNames = T)

addWorksheet(wb, "p.value")
writeData(wb, "p.value", p.mat, rowNames = T)


saveWorkbook(wb, opt$out_file)
