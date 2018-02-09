library(data.table)
if (!"tcga" %in% ls()) {
tcga = as.data.frame(fread("TCGA_normal_primary_zscore.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
}