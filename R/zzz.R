.First.lib <- function (lib, pkg)
{
    library.dynam("tseries", pkg, lib)
    if(!require(ts, quietly = TRUE))
        stop("Package ts is needed.  Stopping")
    mylib <- dirname(.path.package("tseries"))
    ver <- package.description("tseries", lib = mylib)["Version"]
    vertxt <- paste("\n\t`tseries' version:", ver, "\n")
    introtxt <-
        paste("\n\t`tseries' is a package for time series analysis ",
              "with emphasis\n",
              "\ton non-linear modelling.\n",
              "\tSee `library (help=tseries)' for details.\n\n",
              sep = "")
    if(interactive() || .Options$verbose)
        cat(paste (vertxt, introtxt))
}
