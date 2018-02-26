gaiaCNVplot <-  function(calls, cancer=NULL, threshold=0.01)
{
  Calls <- calls[order(calls[,"Region Start[bp]"]),]
  Calls <- Calls[order(Calls[,"Chromosome"]),]
  rownames(Calls) <- NULL
  Chromo <- Calls[,"Chromosome"]
  Gains <-  apply(Calls,1, function(x) ifelse(x["Aberration Kind"]==1, x["score"], 0))
  Losses <- apply(Calls,1,function(x) ifelse(x["Aberration Kind"]==0, x["score"], 0))
  plot(Gains, ylim = c(-max(Calls[,"score"]+2), max(Calls[,"score"]+2)), type  = "h", col  = "red", xlab = "Chromosome", ylab = "Score", xaxt = "n")
  points(-(Losses), type = "h", col = "blue")
  abline(h = 0, cex = 4)
  abline(h = -log10(threshold), col = "orange", cex = 4, main="test")
  abline(h = log10(threshold), col = "orange", cex = 4, main="test")
  uni.chr <- unique(Chromo)
  temp <- rep(0, length(uni.chr))
  for (i in 1 :length(uni.chr)) {
    temp[i] <- max(which(uni.chr[i] == Chromo)) 
  }
  for (i in 1:length(temp)) {
    abline(v = temp[i], col = "black", lty = "dashed", )
  }
  nChroms <- length(uni.chr)
  begin <- c()
  for (d in 1: nChroms) {
    chrom <- sum(Chromo == uni.chr[d])
    begin <- append(begin, chrom)
  }
  temp2 <- rep(0, nChroms)
  for(i in 1:nChroms) {
    if (i == 1) {
      temp2[1] <-(begin[1] * 0.5)
    }
    else  if(i > 1) {
      temp2[i] <- temp[i - 1] + (begin[i] * 0.5)
    }
  }
  uni.chr[uni.chr ==23] <- "X"
  uni.chr[uni.chr ==24] <- "Y"
  for (i in 1:length(temp)) {
    axis(1, at = temp2[i], labels = uni.chr[i], cex.axis = 1) 
  }
  legend(x=1, y=max(Calls[,"score"]+2), y.intersp=0.8, c("Amp" ), pch=15, col=c("red"),  text.font=3)
  legend(x=1,y=-max(Calls[,"score"]+0.5), y.intersp=0.8, c("Del"), pch=15, col=c("blue"), text.font=3)
}

for(cancer in c("SKCM")){
  message(paste0("Starting", cancer))
  # Prepare CNV matrix
  cnvMatrix <- get(load(paste0 (cancer,"_cnv.rda")))
  
  # Add label (0 for loss, 1 for gain)
  cnvMatrix <- cbind(cnvMatrix, Label=NA)
  cnvMatrix[cnvMatrix[,"Segment_Mean"] < -0.3, "Label" ] <- 0
  cnvMatrix[cnvMatrix[,"Segment_Mean"] > 0.3,"Label"] <- 1
  cnvMatrix <- cnvMatrix[!is.na(cnvMatrix$Label),]
  
  # Remove " Segment_Mean" and change col.names
  cnvMatrix <-cnvMatrix[,-6]
  colnames(cnvMatrix) <- c( "Sample.Name", "Chromosome", "Start", "End", "Num.of.Markers", "Aberration")
  
  # Substitute Chromosomes "X" and "Y" with "23" and "24"
  xidx <- which(cnvMatrix$Chromosome=="X")
  yidx <- which(cnvMatrix$Chromosome=="Y")
  cnvMatrix[xidx,"Chromosome"] <- 23
  cnvMatrix[yidx,"Chromosome"] <- 24
  cnvMatrix$Chromosome <- sapply(cnvMatrix$Chromosome,as.integer)
  
  # Recurrent CNV identification with GAIA
  
  # Retrieve probes meta file from broadinstitute website
  # Recurrent CNV identification with GAIA
  gdac.root <- "ftp://ftp.broadinstitute.org/pub/GISTIC2.0/hg19_support/"
  file <- paste0(gdac.root, "genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1.txt")
  # Retrieve probes meta file from broadinstitute website
  if(!file.exists(basename(file))) download(file, basename(file))
  markersMatrix <- readr::read_tsv(basename(file), col_names = FALSE, col_types = "ccn", progress = TRUE)
  colnames(markersMatrix) <- c("Probe.Name", "Chromosome", "Start")
  unique(markersMatrix$Chromosome)
  xidx <- which(markersMatrix$Chromosome=="X")
  yidx <- which(markersMatrix$Chromosome=="Y")
  markersMatrix[xidx,"Chromosome"] <- 23
  markersMatrix[yidx,"Chromosome"] <- 24
  markersMatrix$Chromosome <- sapply(markersMatrix$Chromosome,as.integer)
  markerID <- apply(markersMatrix,1,function(x) paste0(x[2],":",x[3]))
  print(table(duplicated(markerID)))
  print(table(duplicated(markersMatrix$Probe.Name)))
  
  markersMatrix <- markersMatrix[-which(duplicated(markerID)),]
  # Filter markersMatrix for common CNV
  markerID <- apply(markersMatrix,1,function(x) paste0(x[2],":",x[3]))
  
  file <- paste0(gdac.root, "CNV.hg19.bypos.111213.txt")
  if(!file.exists(basename(file))) download(file, basename(file))
  commonCNV <- readr::read_tsv(basename(file), progress = TRUE)
  commonID <- apply(commonCNV,1,function(x) paste0(x[2],":",x[3]))
  print(table(commonID %in% markerID))
  print(table(markerID %in% commonID))
  markersMatrix_fil <- markersMatrix[!markerID %in% commonID,]
  
  markers_obj <- load_markers(as.data.frame(markersMatrix_fil))
  nbsamples <- length(get(paste0("query.",tolower(cancer), ".nocnv"))$results[[1]]$cases)
  cnv_obj <- load_cnv(cnvMatrix, markers_obj, nbsamples)
  results <- runGAIA(cnv_obj, markers_obj, output_file_name=paste0("GAIA_",cancer,"_flt.txt"),
                     aberrations = -1,chromosomes = -1, num_iterations = 10,threshold = 0.25)
  
  # Set q–value threshold
  threshold <- 0.0001
  
  # Plot the results
  RecCNV <- t(apply(results,1,as.numeric))
  colnames(RecCNV) <- colnames(results)
  RecCNV <- cbind(RecCNV, score=0)
  minval <- format(min(RecCNV[RecCNV[,"q–value"]!=0,"q–value"]),scientific =FALSE)
  minval <- substring(minval,1, nchar(minval)-1)
  RecCNV[RecCNV[,"q–value"]==0,"q–value"] <- as.numeric(minval)
  RecCNV[,"score"] <- sapply(RecCNV[,"q–value"],function(x) -log10(as.numeric(x)))
  RecCNV[RecCNV[,"q-value"]==as.numeric(minval),]
  
  gaiaCNVplot(RecCNV,cancer,threshold)
  #save(results, RecCNV, threshold, file = paste0(cancer,"_CNV_results.rda"))
  message(paste0("Results saved as:",cancer,"_CNV_results.rda"))
}