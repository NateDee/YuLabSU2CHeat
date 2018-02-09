library(data.table)
if (!"dbgap" %in% ls()) {
dbgap = as.data.frame(fread("dbgap_su2c_zscore.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
}