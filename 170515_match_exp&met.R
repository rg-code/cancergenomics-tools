i = 0
j = 0
samples.info <- colData(data1)
samples.info$newCol <- NA
for (i in 1:length(testlo1))
{
  for (j in 1:length(data1$"patient"))
  {
    if ((grepl(testlo1[i], data1$"patient"[j])) == TRUE)
    {
      samples.info$newCol[j] <- "DNMT1.low"
  
    }
  }
}
for (i in 1:length(testhigh))
{
  for (j in 1:length(data1$"patient"))
  {
    if ((grepl(testhigh[i], data1$"patient"[j])) == TRUE)
    {
      samples.info$"newCol"[j] <- "DNMT1.high"
    }
  }
}
mid.quantile.samples <- substr(mid.quantile.samples,start = 1, stop = 12)
for (i in 1:length(mid.quantile.samples))
{
  for (j in 1:length(data1$"patient"))
  {
    if ((grepl(mid.quantile.samples[i], data1$"patient"[j])) == TRUE)
    {
      samples.info$"newCol"[j] <- "DNMT1.mid"
    }
  }
}