.First.lib <- function (lib, pkg)
{
  library.dynam("tseries", pkg, lib)
  if (!require (ts, quietly=TRUE))
    stop (paste("Package ts is needed (It should be installed automatically\n",
                "      for version 0.65.0 and later). Stopping"))
  provide(tseries)
  if (interactive() || .Options$verbose)
    cat("\n      `tseries' is a library for time series analysis with emphasize\n",
        "     on non-linear modelling.\n",
        "     See `library (help=tseries)' for details.\n\n")
}
  
