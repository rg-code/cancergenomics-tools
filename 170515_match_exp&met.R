i = 0
j = 0
samples.info <- colData(met)
samples.info$newCol <- NA
low.quantile.samples <- substr(low.quantile.samples,start = 1, stop = 12)
for (i in 1:length(low.quantile.samples))
{
  for (j in 1:length(met$"patient"))
  {
    if ((grepl(low.quantile.samples[i], met$"patient"[j])) == TRUE)
    {
      samples.info$newCol[j] <- "DNMT1.low"
  
    }
  }
}
high.quantile.samples <- substr(high.quantile.samples,start = 1, stop = 12)
for (i in 1:length(high.quantile.samples))
{
  for (j in 1:length(met$"patient"))
  {
    if ((grepl(high.quantile.samples[i], met$"patient"[j])) == TRUE)
    {
      samples.info$"newCol"[j] <- "DNMT1.high"
    }
  }
}
mid.quantile.samples <- substr(mid.quantile.samples,start = 1, stop = 12)
for (i in 1:length(mid.quantile.samples))
{
  for (j in 1:length(met$"patient"))
  {
    if ((grepl(mid.quantile.samples[i], met$"patient"[j])) == TRUE)
    {
      samples.info$"newCol"[j] <- "DNMT1.mid"
    }
  }
}