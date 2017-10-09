
x <- read.csv("G:\\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt",header=TRUE,row.names = 1,sep = "\t")
x <- as.matrix(x)
x[is.nan(x)] <- -5
heatmap.2( x,col =redgreen,
           scale = "none",
           cexRow=0.3, cexCol=0.8,
           margins=c(6,6), trace="none",
           dendrogram='none')






