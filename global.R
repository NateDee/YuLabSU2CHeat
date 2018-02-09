library(data.table)
dbgap = as.data.frame(fread("dbgap_su2c_zscore.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))