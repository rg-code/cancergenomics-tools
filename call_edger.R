source("http://bioconductor.org/biocLite.R")
biocLite("DESeq")
biocLite("edgeR")
library("DESeq")
library("edgeR")
setwd("~/tcga/diffex/edger/")
run_edger = function(dis_state1, dis_state2)
{
  pc <- read.csv("primary_raw_blca.csv",sep=",",header = FALSE)
  nc <- read.csv("normal_raw_blca.csv", sep=",", header = FALSE)
  #tm <- read.csv("metastatic_raw_blca.csv", sep=",",header = FALSE)
  blca <- read.table("edger_blca.csv", sep=",", header = FALSE, row.names = 1)
  #blca <- read.table("edger_blca_nttm.csv", sep=",", header = FALSE, row.names = 1)
  blca_subset <- blca[rowSums(blca)>10,]
  group <- c(rep(dis_state1,ncol(nc)),rep(dis_state2,ncol(pc)))
  #group <- c(rep(dis_state1,ncol(nc)),rep(dis_state2,ncol(tm)))
  d <- DGEList(counts = blca_subset, group=group)
  d <- calcNormFactors(d)
  d <- estimateCommonDisp(d,verbose=T)
  de.tgw <- exactTest(d)
  names(d)
  names(de.tgw)
  head(blca_subset)
  head(de.tgw$table)
  blca_pvalues <- cbind(blca_subset[,1], de.tgw$table)
  blca_pvalues$PValue_fdr <- p.adjust(method="BH", p=blca_pvalues$PValue)
  blca_pvalues_sort <- blca_pvalues[order(blca_pvalues$PValue, blca_pvalues$PValue_fdr),]
  #View(blca_pvalues_sort)
  
  #collect data
  blca_genes_in_data_set <- rownames(blca_pvalues_sort[0])
  blca_DE_0.0001_genes <- blca_pvalues_sort[blca_pvalues_sort$PValue<0.0001]
  blca_DE_0.0001_genes_names <- rownames(blca_DE_0.001_genes[0])
  
  #write stuff
  write.csv(blca_DE_0.001_genes,"blca_DE_0.0001_genes.csv")
  #write.csv(blca_DE_0.001_genes,"blca_DE_0.001_genes_nttm.csv")
  #write.csv(blca_DE_0.001_genes_names, "blca_DE_0.001_genes_names_nttm.csv")
  write.csv(blca_DE_0.001_genes_names, "blca_DE_0.0001_genes_names.csv")
}
