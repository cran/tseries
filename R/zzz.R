.onAttach <-
function(lib, pkg)
{
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
