# Copyright (C) 1997-1999  Adrian Trapletti
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, write to the Free
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#
# tsparam object
#


tsparam <- function (x, y, xlab="", ylab="", lab="", ybnd=NULL, ysgnf=NULL)
{
  if (!is.vector(x) | !is.vector(y))
    stop ("x and y must be vectors")
  if (length(x) != length(y))
    stop ("x and y have not the same length")
  if (!is.character(xlab) | !is.character(ylab) | !is.character(lab))
    stop ("xlab, ylab, and lab must be strings")
  if ((!is.null(ybnd)) & (!is.list(ybnd)))
    stop ("ybnd must be a list")
  if ((!is.null(ybnd)) & (any(sapply(ybnd,length) != length(x))))
    stop ("x and the elements of ybnd must have the same length")
  if ((!is.null(ysgnf)) & !is.character(ysgnf))
    stop ("ysgnf must be a vector of strings")
  if ((!is.null(ysgnf)) & (length(x) != length(ysgnf)))
    stop ("x and ysgnf must have the same length")
  structure (list(x=x, y=y), xlab=xlab, ylab=ylab, ybnd=ybnd,
             ysgnf=ysgnf, lab=lab, class="tsparam")
}

plot.tsparam <- function (obj, xlim = range(obj$x), ylim = range(attr(obj,"ybnd"),obj$y),
                          xlab = attr(obj,"xlab"), ylab = attr(obj,"ylab"),
                          type = "l", main = "", grid = TRUE, col = "black", 
                          bcol = "blue", gridcol = "lightgrey") 
{
  if (!inherits(obj, "tsparam")) 
    stop ("obj is not of class tsparam")
  if (is.complex(obj$x) | is.complex(obj$y))
    stop ("complex valued obj can not be plotted")
  plot (obj$x, obj$y, xlim = xlim, ylim = ylim, pch = 0, cex = 0.5,
        main = main, xlab = xlab, ylab = ylab, type = type, col = col)
  if (attr(obj,"lab") != "")
    mtext(paste("[",attr(obj,"lab"),"]"), side = 1, line = 4, cex = 0.8)
  for (i in seq (along = attr(obj, "ybnd")))
    lines (obj$x, attr(obj, "ybnd")[[i]], lty = 2, col = bcol)
  if (grid)
  {
    xaxp <- par("xaxp")
    yaxp <- par("yaxp")
    xseq <- seq(xaxp[1], xaxp[2], length = xaxp[3] + 1)
    yseq <- seq(yaxp[1], yaxp[2], length = yaxp[3] + 1)
    abline(v = xseq, lty = 3, col = gridcol)
    abline(h = yseq, lty = 3, col = gridcol)
    plot.xy(xy.coords(obj$x, obj$y), pch = 0, cex = 0.5, 
            type = type, col = col)
  }
  invisible(obj)
}

print.tsparam <- function (obj, digits = max(3,.Options$digits-3), ...)
{
  if (!inherits(obj, "tsparam")) 
    stop ("obj is not of class tsparam")
  if (is.null(attr(obj,"ysgnf")))
  {
    pr <- cbind (format(obj$x, digits=digits),
                 format(obj$y, digits=digits))
    colnames(pr) <- c (abbreviate(attr(obj,"xlab"),6),
                       abbreviate(attr(obj,"ylab"),6))
  }
  else
  {
    pr <- cbind (format(obj$x, digits=digits),
                 format(obj$y, digits=digits),
                 formatC(attr(obj,"ysgnf"), width=3))
    colnames(pr) <- c (abbreviate(attr(obj,"xlab"),6),
                       abbreviate(attr(obj,"ylab"),6),
                       abbreviate("significance",6))
  }
  rownames(pr) <- rep ("", length(obj$x))
  print (pr, quote = F, ...)
  if (attr(obj,"lab") != "")
  {
    cat ("---\n")
    cat (attr(obj,"lab"),"\n")
  }
  invisible (obj)
}







