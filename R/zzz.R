.First.lib <-
function(lib, pkg)
{
    library.dynam("tseries", pkg, lib)
    require("stats", quietly = TRUE)
    mylib <- dirname(system.file(package = "tseries"))
    ver <- packageDescription("tseries", lib = mylib)["Version"]
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
