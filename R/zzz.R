.First.lib <- function (lib, pkg)
{
  library.dynam("tseries", pkg, lib)
  if (!require (ts, quietly=TRUE))
    stop ("Package ts is needed. Stopping")
  filenm <- paste (.lib.loc, "/tseries/DESCRIPTION", sep="")
  descidx <- file.exists (filenm)
  descfn <- paste (filenm)[descidx]
  desc <- scan (descfn, what=character(), quiet = TRUE)
  ver <- desc[which (desc=="Version:")+1]
  vertxt <- paste ("\n      `tseries' version:", ver, "\n")
  introtxt <- paste ("\n      `tseries' is a package for time series analysis with emphasize\n",
                     "       on non-linear modelling.\n",
                     "       See `library (help=tseries)' for details.\n\n", sep="")
  if (interactive() || .Options$verbose)
    cat (paste (vertxt, introtxt))
}
  
