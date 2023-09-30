# Jaccard index with p value
# Implemented from Kumar et al. PLoS One (2017) doi: 10.1371/journal.pone.0187132
# P value is computed as one sided

library(optparse)
library(tidyverse)
library(openxlsx)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="OTU table (.csv).", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", 
              help="Output file for the Jaccard similarity with the p value.", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

jaccard <- function(x, y){
  xy <- x + y
  x_intersect_y <- sum(xy == 2)
  x_union_y <- sum(xy > 0)
  
  j.stat <- x_intersect_y / x_union_y
  
  m.value <- qhyper(0.5, sum(x), (length(x)-sum(x)), sum(y))
  
  if(x_intersect_y <= m.value){
    sign.value <- (-1)
    p.value <- phyper(x_intersect_y, sum(x), (length(x)-sum(x)), sum(y), lower.tail = T)
    j.stat <- (1-j.stat)*(-1)
  } else {
    sign.value <- (+1)
    p.value <- phyper(x_intersect_y, sum(x), (length(x)-sum(x)), sum(y), lower.tail = F)
  }
  
  return(c(j.stat, p.value, sign.value))
}

otu.table <- read.csv(opt$file, row.names = 1)

## Make the symmetrical matrices that will store the results (#OTUs x #OTUs)

J.mat <- matrix(data = 0,
                nrow = nrow(otu.table),
                ncol = nrow(otu.table)) %>%
  `colnames<-`(rownames(otu.table)) %>%
  `rownames<-`(rownames(otu.table))

p.mat <- matrix(data = 0,
                nrow = nrow(otu.table),
                ncol = nrow(otu.table)) %>%
  `colnames<-`(rownames(otu.table)) %>%
  `rownames<-`(rownames(otu.table))

#s.mat <- matrix(data = 0,
#                nrow = nrow(otu.table),
#                ncol = nrow(otu.table)) %>%
#  `colnames<-`(rownames(otu.table)) %>%
#  `rownames<-`(rownames(otu.table))

## Compute Jaccard similarity, p value and sign for the upper triangular matrices

N <- ((nrow(otu.table)**2)-nrow(otu.table))/2
print(paste0("Computing ", N, " jaccard similarities"))

counter <- 0
for(i in 1:(nrow(otu.table)-1)){
  for(j in (i+1):nrow(otu.table)){
    J.val <- jaccard(otu.table[i, ], otu.table[j, ])
    #print(c(rownames(otu.table)[i], rownames(otu.table)[j]))
    #print(c(sum(is.na(otu.table[i, ])), sum(is.na(otu.table[j, ]))))
    #print(J.val)
    J.mat[i, j] <- J.val[1]
    p.mat[i, j] <- J.val[2]
    #s.mat[i, j] <- J.val[3]
    
    counter <- counter + 1
    if((counter%%100000)==0){
      print(paste0("Jaccard similarity ", counter, " out of ", N, ": ", round((counter/N)*100, digits = 2), "% done."))
    }
    
  }
}

## Since the similarity measure is symmetrical add the transposed matrix

J.mat <- J.mat + t(J.mat)
p.mat <- p.mat + t(p.mat)
#s.mat <- s.mat + t(s.mat)

## Set the diagonal for the identity

diag(J.mat) <- 1
diag(p.mat) <- 0
#diag(s.mat) <- +1

if(file.exists(opt$out_file)){
  file.remove(opt$out_file)
}

wb <- createWorkbook()

addWorksheet(wb, "jaccard")
writeData(wb, "jaccard", J.mat, rowNames = T)

addWorksheet(wb, "p.value")
writeData(wb, "p.value", p.mat, rowNames = T)

#addWorksheet(wb, "sign")
#writeData(wb, "sign", s.mat, rowNames = T)

saveWorkbook(wb, opt$out_file)
