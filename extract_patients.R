files <- list.files("~/tcga/mrnaseq_all/all_raw/")
setwd("~/tcga/mrnaseq_all/all_raw/")
#create normal list
normal <- file("~/tcga/mrnaseq_all/normal_mrnaseq.txt")
mylist_n <- strsplit(readLines(normal), "\n")
close(normal)

#create tp list
tp <- file("~/tcga/mrnaseq_all/primary_solid_mrnaseq.txt")
mylist_tp <- strsplit(readLines(tp), "\n")
close(tp)

#create tm list
tm <- file("~/tcga/mrnaseq_all/metastatic_mrnaseq.txt")
mylist_tm <- strsplit(readLines(tm), "\n")
close(tm)

met_file_raw = c()
primary_file_raw = c()
normal_file_raw = c()

#loop for all the cancer types
for (i in files) {
  test <- read.delim(i, header=T)
  colnum <- ncol(test)
  colnam <- colnames(test)
  
  #classify patients across all cancer
  for ( j in colnam ) {
    #metastatic patients
    for ( k in mylist_tm ) {
      if ((grepl(j,k)) == TRUE ) {
        extract_tm <- (test[,j])
        extract_tm <- as.matrix(extract_tm)
        met_file_raw = c(met_file_raw, extract_tm)
      }
    }
    
    #normal patients
    for ( m in mylist_n ) {
      if ((grepl(j,m)) == TRUE ) {
        extract_n <- (test[,j])
        extract_n <- as.matrix(extract_n)
        normal_file_raw = c(normal_file_raw, extract_n)
      }
    }
    
    #primary solid patients
    for ( n in mylist_tp ) {
      if ((grepl(j,n)) == TRUE ) {
        extract_tp <- (test[,j])
        extract_tp <- as.matrix(extract_tp)
        primary_file_raw <- c(primary_file_raw, extract_tp)
      }
    }
  }
}

# write matrix to file
setwd("~/tcga/mrnaseq_all/all_raw/")
WriteMatrix = function(data, fnm, nrows = NA, ncols = NA)
{
  if (is.na(nrows))
  {
    nrows = nrow(data)
  }
  if (is.na(ncols))
  {
    ncols = ncol(data)
  }
  sink(file = fnm)
  for (i in 1:nrows)
  {
    cat(sprintf("%8.3f", data[i]))
    for (j in 2:ncols)
    {
      cat(sprintf(", %8.3f", data[(j - 1) * nrows + i]))
    }
    cat("\n")
  }
  sink()
}
