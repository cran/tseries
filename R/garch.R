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
# GARCH class
#


garch <- function (x, order = c(1, 1), coef = NULL, itmax = 200, eps = NULL,
                   grad = c("analytical","numerical"), series = NULL, trace = TRUE, ...)
{
  if (NCOL(x) > 1) stop ("x is not a vector or univariate time series")
  if (!is.vector(order)) stop ("order is not a vector")
  grad <- match.arg (grad)
  switch (grad,
          analytical = (agrad <- TRUE),
          numerical = (agrad <- FALSE))
  if (is.null(series)) series <- deparse(substitute(x))
  ists <- is.ts(x)
  x <- as.ts(x)
  xfreq <- frequency(x)
  if (any(is.na(x))) stop ("NAs in x")
  if (ists) xtsp <- tsp(x)
  x <- as.matrix(x)
  n <- nrow(x)
  e <- double(n)
  ncoef <- order[1]+order[2]+1
  hess <- matrix (0.0, ncoef, ncoef)
  small <- 0.05
  if (is.null(coef)) coef <- c(var(x)*(1.0-small*(ncoef-1)),rep(small,ncoef-1))
  if (!is.vector(coef)) stop ("coef is not a vector")
  if (ncoef != length(coef)) stop ("incorrect length of coef")
  if (is.null(eps)) eps <- Machine()$double.eps
  nlikeli <- 1.0e+10
  fit <- .C ("fit_garch", as.vector(x,mode="double"), as.integer(n),
             coef=as.vector(coef,mode="double"), as.integer(order[1]),
             as.integer(order[2]), as.integer(itmax), as.double(eps),
             nlikeli=as.double(nlikeli), as.integer(agrad), as.integer(trace),
             PACKAGE="tseries")
  pred <- .C ("pred_garch", as.vector(x,mode="double"), e=as.vector(e,mode="double"),
              as.integer(n), as.vector(fit$coef,mode="double"), as.integer(order[1]),
              as.integer(order[2]), PACKAGE="tseries")
  com.hess <- .C ("ophess_garch", as.vector(x,mode="double"), as.integer(n),
                  as.vector(fit$coef,mode="double"), hess=as.matrix(hess),
                  as.integer(order[1]), as.integer(order[2]), PACKAGE="tseries")
  rank <- qr(com.hess$hess,...)$rank
  if (rank != ncoef)
  {
    se.garch <- rep (NA, ncoef)
    cat ("Warning: singular information\n")
  }
  else
    se.garch <- sqrt(diag(solve(com.hess$hess)))
  sigt <- sqrt(pred$e)
  sigt[1:max(order[1],order[2])] <- rep (NA, max(order[1],order[2]))
  f <- cbind(sigt,-sigt)
  colnames(f) <- c("sigt","-sigt")
  e <- as.vector(x)/sigt  
  if (ists)
  {
    attr(e, "tsp") <-  attr(f, "tsp") <- xtsp
    attr(e, "class") <- attr(f, "class") <- "ts"
  }
  names(order) <- c("p","q")
  coef <- fit$coef
  nam.coef <- "a0"
  if (order[2] > 0) nam.coef <- c(nam.coef, paste("a", seq(order[2]), sep = ""))
  if (order[1] > 0) nam.coef <- c(nam.coef, paste("b", seq(order[1]), sep = ""))
  names(coef) <- nam.coef
  names(se.garch) <- nam.coef
  garch <- list (order = order, coef = coef, n.likeli = fit$nlikeli,
                 n.used = n, residuals = e, fitted.values = f, series = series,
                 frequency = xfreq, call = match.call(), asy.se.coef = se.garch)
  class(garch) <- "garch"
  return (garch)
}

coef.garch <- function (object)
{
  if (!inherits(object, "garch")) stop ("method is only for garch objects")
  return (object$coef)
}

residuals.garch <- function (object)
{
  if (!inherits(object, "garch")) stop ("method is only for garch objects")
  return (object$residuals)
}

fitted.garch <- function (object)
{
  if (!inherits(object, "garch")) stop ("method is only for garch objects")
  return (object$fitted.values)
}

print.garch <- function (object, digits = max(3,.Options$digits-3))
{
  if (!inherits(object, "garch")) stop ("method is only for garch objects")
  cat ("\nCall:\n", deparse(object$call), "\n\n", sep = "")
  cat ("Coefficient(s):\n")
  print.default (format(coef(object), digits = digits), print.gap = 2, quote = FALSE)
  cat ("\n")
  invisible (object)
}

summary.garch <- function (object)
{
  if (!inherits(object, "garch")) stop ("method is only for garch objects")
  ans <- NULL
  ans$residuals <- na.remove(object$residuals)
  tval <- object$coef/object$asy.se.coef
  ans$coef <- cbind (object$coef, object$asy.se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(ans$coef) <- list(names(object$coef),
                             c(" Estimate"," Std. Error"," t value","Pr(>|t|)"))
  ans$call <- object$call
  ans$order <- object$order
  Residuals <- ans$residuals
  ans$j.b.test <- jarque.bera.test(Residuals)
  Squared.Residuals <- ans$residuals^2
  ans$l.b.test <- Box.test (Squared.Residuals, type = "Ljung-Box")
  class(ans) <- "summary.garch"
  return (ans)
}

plot.garch <- function (object, ask = interactive())
{
  if (!inherits(object, "garch")) stop ("method is only for garch objects")
  op <- par()
  par (ask = ask, mfrow=c(2,1))
  x <- eval (parse(text=object$series))
  if (any(is.na(x))) stop ("NAs in x")
  plot(x, main = object$series, ylab = "Series")
  plot(object$residuals, main = "Residuals", ylab = "Series")
  hist (x, main = paste("Histogram of",object$series), xlab = "Series")
  hist (object$residuals, main = "Histogram of Residuals", xlab = "Series")
  qqnorm (x, main = paste("Q-Q Plot of",object$series), xlab = "Normal Quantiles")
  qqnorm (object$residuals, main = "Q-Q Plot of Residuals", xlab = "Normal Quantiles")
  acf (x^2, main = paste("ACF of Squared",object$series))
  acf (object$residuals^2, main = "ACF of Squared Residuals", na.action=na.remove)
  par (ask = op$ask, mfrow = op$mfrow)
  invisible (object)
}

print.summary.garch <- function (object, digits = max(3,.Options$digits-3),
                                 signif.stars = .Options$show.signif.stars, ...)
{
  if (!inherits(object, "summary.garch")) stop ("method is only for summary.garch objects")
  cat ("\nCall:\n", deparse(object$call), "\n", sep = "")
  cat ("\nModel:\nGARCH(", object$order[1], ",", object$order[2], ")", "\n", sep = "")
  cat ("\nResiduals:\n")
  rq <- structure(quantile(object$residuals), names = c("Min","1Q","Median","3Q","Max"))
  print (rq, digits=digits, ...)
  cat ("\nCoefficient(s):\n")
  print.coefmat (object$coef, digits = digits, signif.stars = signif.stars, ...)
  cat ("\nDiagnostic Tests:")
  print (object$j.b.test)
  print (object$l.b.test)
  invisible (object)
}

predict.garch <- function (object, newdata, genuine = FALSE)
{
  if (!inherits(object, "garch")) stop ("method is only for garch objects")
  if (missing(newdata))
  {
    newdata <- eval (parse(text=object$series))
    if (any(is.na(newdata))) stop ("NAs in newdata")
  }
  if (NCOL(newdata) > 1) stop ("newdata is not a vector or univariate time series")
  ists <- is.ts(newdata)
  if (ists) newdata.tsp <- tsp(newdata)
  newdata <- as.matrix(newdata)
  n <- nrow(newdata)
  if (genuine) h <- double(n+1)
  else h <- double(n)
  pred <- .C ("pred_garch", as.vector(newdata,mode="double"), h=as.vector(h,mode="double"),
              as.integer(n), as.vector(object$coef,mode="double"), as.integer(object$order[1]),
              as.integer(object$order[2]), PACKAGE="tseries")
  pred$h <- sqrt(pred$h)
  pred$h[1:max(object$order[1],object$order[2])] <- rep (NA, max(object$order[1],object$order[2]))
  pred$h <- cbind(pred$h,-pred$h)
  if (ists)
  {
    if (genuine) attr(pred$h, "tsp") <- c(newdata.tsp[1],
                                          newdata.tsp[2]+1/newdata.tsp[3],
                                          newdata.tsp[3])
    else attr(pred$h, "tsp") <- newdata.tsp
    attr(pred$h, "class") <- "ts"
  }
  return (pred$h)
}

