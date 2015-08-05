# data:  Can be a matrix or vector.
# fnm:   Output file name.
# nrows: Number of rows in output matrix (lines).
# ncols: Number of cols in output matrix.
#
# Precision can be altered by changing the %8.3f.
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