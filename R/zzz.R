.First.lib <-
function(lib, pkg)
{
    library.dynam("tseries", pkg, lib)
    RisLT19 <- ((R.version$major == 1)
                && (as.numeric(R.version$minor) < 9))
    if(RisLT19) {
        if(!require("ts", quietly = TRUE))
            stop("Package", sQuote("ts"), "is needed.  Stopping")
    } else {
        require("stats", quietly = TRUE)
    }
    mylib <- dirname(system.file(package = "tseries"))
    ver <- if(RisLT19)
        package.description("tseries", lib = mylib)["Version"]
    else
        packageDescription("tseries", lib = mylib)["Version"]
    txt <- c("\n",
             paste(sQuote("tseries"), "version:", ver),
             "\n",
             paste(sQuote("tseries"),
                   "is a package for time series analysis",
                   "and computational finance."),
             "\n",
             paste("See",
                   sQuote("library(help=\"tseries\")"),
                   "for details."),
             "\n")
    if(interactive() || getOption("verbose"))
        writeLines(strwrap(txt, indent = 4, exdent = 4))
}
