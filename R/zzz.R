.First.lib <- function (lib, pkg)
{
  library.dynam("tseries", pkg, lib)
  if (!require (ts, quietly=TRUE))
    stop ("Package ts is needed. Stopping")
  mylib <- .path.package("tseries")
  mylib <- substr(mylib, 1, nchar(mylib)-7)
  ver <- package.description("tseries", lib=mylib)$Version
  vertxt <- paste ("\n      `tseries' version:", ver, "\n")
  introtxt <- paste ("\n      `tseries' is a package for time series analysis with emphasize\n",
                     "       on non-linear modelling.\n",
                     "       See `library (help=tseries)' for details.\n\n", sep="")
  if (interactive() || .Options$verbose)
    cat (paste (vertxt, introtxt))
}
  
