# Copyright (C) 1997-2000  Adrian Trapletti
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
# ARMA class
#


arma <- function (x, order = c(1, 1), lag = NULL, coef = NULL, include.intercept = TRUE,
                  series = NULL, qr.tol = 1e-07, ...)
{
  
  seqN <- function(N) {
    if (0==length(N)) NULL else if (N<=0) NULL else seq(N)
  }
  
  err <- function (coef)
  {
    u <- double(n)
    u[seqN(max.order)] <- 0
    u <- .C("arma", as.vector(x, mode="double"), u=as.vector(u),
            as.vector(coef, mode="double"), as.integer(lag$ar), as.integer(lag$ma),
            as.integer(ar.l), as.integer(ma.l), as.integer(max.order), as.integer(n),
            as.integer(include.intercept), PACKAGE="tseries")$u
    return (sum(u^2))
  }
  
  resid <- function (coef)
  {
    u <- double(n)
    u[seqN(max.order)] <- 0
    u <- .C("arma", as.vector(x, mode="double"), u=as.vector(u),
            as.vector(coef, mode="double"), as.integer(lag$ar), as.integer(lag$ma),
            as.integer(ar.l), as.integer(ma.l), as.integer(max.order), as.integer(n),
            as.integer(include.intercept), PACKAGE="tseries")$u
    return (u)
  }
  
  arma.init <- function ()
  {
    k <- round(1.1*log(n))
    e <- na.omit (drop (ar.ols (x, order.max=k, aic=FALSE, demean=FALSE,
                                intercept=include.intercept)$resid))
    ee <- embed (e, max.order+1)
    xx <- embed (x[-(1:k)], max.order+1)
    if (include.intercept == TRUE)
    {
      if (is.null(lag$ar))
        coef <- lm (xx[,1]~ee[,lag$ma+1])$coef
      else if (is.null(lag$ma))
        coef <- lm (xx[,1]~xx[,lag$ar+1])$coef
      else
        coef <- lm (xx[,1]~xx[,lag$ar+1]+ee[,lag$ma+1])$coef
      coef <- c (coef[-1], coef[1])
    } 
    else 
    {
      if (is.null(lag$ar))
        coef <- lm (xx[,1]~ee[,lag$ma+1]-1)$coef
      else if (is.null(lag$ma))
        coef <- lm (xx[,1]~xx[,lag$ar+1]-1)$coef
      else
        coef <- lm (xx[,1]~xx[,lag$ar+1]+ee[,lag$ma+1]-1)$coef
    }
    return (coef) 
  }
  
  if (!is.null (order) & !is.null (lag))
    warning ("order is ignored")
  if (is.null (order) & is.null (lag))
    stop ("order or lag must be given")
  if (is.null(lag) & !is.null(order))
    lag <- list (ar=seqN(order[1]), ma=seqN(order[2]))
  lag$ar <- unique(lag$ar)
  lag$ma <- unique(lag$ma)
  max.order <- max(unlist(lag),0)
  ar.l <- length(lag$ar)
  ma.l <- length(lag$ma)
  if (NCOL(x) > 1) stop ("x is not a vector or univariate time series")
  if (is.null(series)) series <- deparse(substitute(x))
  ists <- is.ts(x)
  x <- as.ts(x)
  xfreq <- frequency(x)
  if (any(is.na(x))) stop ("NAs in x")
  if (ists) xtsp <- tsp(x)
  n <- length(x)
  if (!is.null(unlist(lag)))
    if (min(unlist(lag)) < 1 | max(unlist(lag)) > (n-1))
      stop ("invalid lag")
  ncoef <- length(unlist(lag))+as.numeric(include.intercept)
  if (is.null(coef))
  {
    if (!is.null(unlist(lag)))
      coef <- arma.init ()
    else
      coef <- 0
  }
  if (length(coef) != ncoef) stop ("invalid coef")
  md <- optim (coef, err, gr=NULL, hessian=TRUE, ...)
  coef <- md$par
  rank <- qr(md$hessian, qr.tol)$rank
  if (rank != ncoef)
  {
    se <- rep (NA, ncoef)
    cat ("Warning: singular Hessian\n")
  }
  else
  {
    di <- diag (2*md$value/n*solve(md$hessian))
    if (any (di < 0)) 
      cat ("Warning: Hessian negative-semidefinite\n")
    se <- sqrt (di)
  }
  e <- resid (coef)
  e[seqN(max.order)] <- NA
  f <- x-e
  if (ists)
  {
    attr(e, "tsp") <- xtsp
    attr(e, "class") <- "ts"
    attr(f, "tsp") <- xtsp
    attr(f, "class") <- "ts"
  }
  if (!is.null(lag$ar)) nam.ar <- paste("ar", lag$ar, sep = "")
  else nam.ar <- NULL
  if (!is.null(lag$ma)) nam.ma <- paste("ma", lag$ma, sep = "")
  else nam.ma <- NULL
  if (include.intercept) nam.int <- "intercept"
  else nam.int <- NULL
  nam.coef <- c(nam.ar, nam.ma, nam.int)
  names(coef) <- nam.coef
  names(se) <- nam.coef
  arma <- list (coef = coef, css = md$value, n.used = n,
                residuals = e, fitted.values = f, series = series, frequency = xfreq,
                call = match.call(), asy.se.coef = se, lag = lag,
                convergence = md$convergence, include.intercept = include.intercept)
  class(arma) <- "arma"
  return(arma)
}

coef.arma <- function (object)
{
  if (!inherits(object, "arma")) stop ("method is only for arma objects")
  return (object$coef)
}

residuals.arma <- function (object)
{
  if (!inherits(object, "arma")) stop ("method is only for arma objects")
  return (object$residuals)
}

fitted.arma <- function (object)
{
  if (!inherits(object, "arma")) stop ("method is only for arma objects")
  return (object$fitted.values)
}

print.arma <- function (object, digits = max(3,.Options$digits-3))
{
  if (!inherits(object, "arma")) stop ("method is only for arma objects")
  cat ("\nCall:\n", deparse(object$call), "\n\n", sep = "")
  cat ("Coefficient(s):\n")
  print.default (format(coef(object), digits = digits), print.gap = 2, quote = FALSE)
  cat ("\n")
  invisible (object)
}

summary.arma <- function (object)
{
  if (!inherits(object, "arma")) stop ("method is only for arma objects")
  ans <- NULL
  ans$residuals <- na.remove(object$residuals)
  tval <- object$coef/object$asy.se.coef
  ans$coef <- cbind (object$coef, object$asy.se.coef, tval, 2*(1-pnorm(abs(tval))))
  dimnames(ans$coef) <- list(names(object$coef),
                             c(" Estimate"," Std. Error"," t value","Pr(>|t|)"))
  ans$call <- object$call
  ans$nn <- object$nn
  ans$css <- object$css
  ans$var <- var(ans$residuals)
  ans$aic <- object$n.used*(1+log(2*pi))+object$n.used*log(ans$var)+2*length(object$coef)
  ans$p <- max(object$lag$ar,0)
  ans$q <- max(object$lag$ma,0)
  class(ans) <- "summary.arma"
  return (ans)
}

print.summary.arma <- function (object, digits = max(3,.Options$digits-3),
                                signif.stars = .Options$show.signif.stars, ...)
{
  if (!inherits(object, "summary.arma")) stop ("method is only for summary.arma objects")
  cat ("\nCall:\n", deparse(object$call), "\n", sep = "")
  cat ("\nModel:\nARMA(",object$p,",",object$q,")\n", sep = "")
  cat ("\nResiduals:\n")
  rq <- structure(quantile(object$residuals), names = c("Min","1Q","Median","3Q","Max"))
  print (rq, digits=digits, ...)
  cat ("\nCoefficient(s):\n")
  print.coefmat (object$coef, digits = digits, signif.stars = signif.stars, ...)
  cat("\nFit:\n")
  cat ("sigma^2 estimated as ", format(object$var, digits = digits), 
       ",  Conditional Sum-of-Squares = ", format(round(object$css, 2)), 
       ",  AIC = ", format(round(object$aic, 2)), "\n", sep = "")
  cat ("\n")
  invisible (object)
}

plot.arma <- function (object, ask = interactive())
{
  if (!inherits(object, "arma")) stop ("method is only for arma objects")
  op <- par()
  par (ask = ask, mfrow=c(2,1))
  x <- eval (parse(text=object$series))
  if (any(is.na(x))) stop ("NAs in x")
  plot(x, main = object$series, ylab = "Series")
  plot(object$residuals, main = "Residuals", ylab = "Series")
  acf (x, main = paste("ACF of",object$series))
  acf (object$residuals, main = "ACF of Residuals", na.action=na.remove)
  pacf (x, main = paste("PACF of",object$series))
  pacf (object$residuals, main = "PACF of Residuals", na.action=na.remove)
  par (ask = op$ask, mfrow = op$mfrow)
  invisible (object)
}

