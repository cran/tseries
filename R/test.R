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
# Time series tests
#


box.test <- function (x, lag = 1, pierce = TRUE)
{
  if (!is.vector(x) & !is.univariate.ts(x))
    stop ("x is not a vector or univariate time series")
  cor <- acf (x, lag = lag, pl = FALSE)
  n <- length(x)
  DNAME <- deparse(substitute(x))
  PARAMETER <- lag
  obs <- cor$y[2:(lag+1)]
  if (pierce)
  {
    METHOD <- "Box-Pierce test"
    STATISTIC <- n*sum(obs^2)
    PVAL <- 1-pchisq(STATISTIC,lag)
  }
  else
  {
    METHOD <- "Box-Ljung test"
    STATISTIC <- n*(n+2)*sum(1/seq(n-1,n-lag)*obs^2)
    PVAL <- 1-pchisq(STATISTIC,lag)
  } 
  names(STATISTIC) <- "X-squared"
  names(PARAMETER) <- "df"
  structure(list(statistic = STATISTIC,
                 parameter = PARAMETER,
		 p.value = PVAL,
		 method = METHOD,
		 data.name = DNAME),
	    class = "htest")
}

runs.test <- function (x)
{
  if (!is.vector(x) & !is.univariate.ts(x))
    stop ("x is not a vector or univariate time series")
  DNAME <- deparse (substitute(x))
  if (any (x == 0.0))
  {
    cat ("Removed", length (x[x==0.0]), "zero(es)\n")
    x <- x[x!=0.0]
  }
  d <- diff (sign(x))
  f <- factor (d)
  sp <- split (d, f)
  resS <- sapply (sp, length)
  resL <- lapply (sp, length)
  n <- length (x)
  sum2 <- sum (resS^2)
  sum3 <- sum (resS^3)
  m <- (n*(n+1)-sum2)/n
  s <- (sum2*(sum2+n*(n+1))-2*n*sum3-n^3)/(n^2*(n-1))
  R <- 1
  if (!is.null(resL$"-2"))
    R <- R+resL$"-2"
  if (!is.null(resL$"2"))
    R <- R+resL$"2"
  STATISTIC <- ((R+0.5)-m)/s
  METHOD <- "Runs test"
  PVAL <- 2 * pnorm (-abs(STATISTIC))
  names(STATISTIC) <- "Standard Normal"
  structure(list(statistic = STATISTIC,
		 p.value = PVAL,
		 method = METHOD,
		 data.name = DNAME),
	    class = "htest")
}

pp.test <- function (x, lshort = TRUE)
{
  if (!is.vector(x) & !is.univariate.ts(x))
    stop ("x is not a vector or univariate time series")
  DNAME <- deparse(substitute(x))
  z <- embed (x, 2)
  yt <- z[,1]
  yt1 <- z[,2]
  n <- length (yt)
  tt <- (1:n)-n/2
  res <- lm (yt~1+tt+yt1)
  res.sum <- summary (res)
  tstat <- (res.sum$coefficients[3,1]-1)/res.sum$coefficients[3,2]
  u <- residuals (res)
  ssqru <- sum(u^2)/n
  if (lshort)
    l <- trunc(4*(n/100)^0.25)
  else
    l <- trunc(12*(n/100)^0.25)
  ssqrtl <- .C ("R_pp_sum", as.vector(u,mode="double"), as.integer(n),
                as.integer(l), trm=as.double(ssqru))
  ssqrtl <- ssqrtl$trm
  n2 <- n^2
  trm1 <- n2*(n2-1)*sum(yt1^2)/12
  trm2 <- n*sum(yt1*(1:n))^2
  trm3 <- n*(n+1)*sum(yt1*(1:n))*sum(yt1)
  trm4 <- (n*(n+1)*(2*n+1)*sum(yt1)^2)/6
  Dx <- trm1-trm2+trm3-trm4
  STAT <- sqrt(ssqru)/sqrt(ssqrtl)*tstat-(n^3)/(4*sqrt(3)*sqrt(Dx)*sqrt(ssqrtl))*(ssqrtl-ssqru)
  table <- cbind(c(4.38,4.15,4.04,3.99,3.98,3.96),
                 c(3.95,3.80,3.73,3.69,3.68,3.66),
                 c(3.60,3.50,3.45,3.43,3.42,3.41),
                 c(3.24,3.18,3.15,3.13,3.13,3.12),
                 c(1.14,1.19,1.22,1.23,1.24,1.25),
                 c(0.80,0.87,0.90,0.92,0.93,0.94),
                 c(0.50,0.58,0.62,0.64,0.65,0.66),
                 c(0.15,0.24,0.28,0.31,0.32,0.33))
  table <- -table
  tablen <- dim(table)[2]
  tableT <- c(25,50,100,250,500,100000)
  tablep <- c(0.01,0.025,0.05,0.10,0.90,0.95,0.975,0.99)
  tableipl <- numeric(tablen)
  for (i in (1:tablen))
    tableipl[i] <- approx (tableT,table[,i],n,rule=2)$y
  PVAL <- approx (tableipl,tablep,STAT,rule=2)$y
  PARAMETER <- l
  METHOD <- "Phillips-Perron Unit Root Test"
  names(STAT) <- "Dickey-Fuller"
  names(PARAMETER) <- "Truncation lag parameter"
  structure(list(statistic = STAT, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name = DNAME), 
            class = "htest")
}

bds.test <- function (x, m = 2, eps = seq(0.5*sd(x),2*sd(x),length=4), trace = FALSE)
{
  if (!is.vector(x) & !is.univariate.ts(x))
    stop ("x is not a vector or univariate time series")
  if (m < 2) stop ("m is less than 2")
  if (any(eps<=0)) stop ("invalid eps")
  DNAME <- deparse(substitute(x))
  n <- length(x)
  k <- length(eps)
  cc <- double(m+1)
  cstan <- double(m+1)
  STATISTIC <- matrix(0,m-1,k)
  for (i in (1:k))
  {
    res <- .C("bdstest_main", as.integer(n), as.integer(m), as.vector(x,mode="double"),
              as.vector(cc), cstan=as.vector(cstan), as.double(eps[i]), as.integer(trace))
    STATISTIC[,i] <- res$cstan[2:m+1]
  }
  colnames(STATISTIC) <- eps
  rownames(STATISTIC) <- 2:m
  PVAL <- 2 * pnorm (-abs(STATISTIC))
  colnames(PVAL) <- eps
  rownames(PVAL) <- 2:m
  METHOD <- "BDS Test"
  PARAMETER <- list (m = 2:m, eps = eps)
  structure(list(statistic = STATISTIC, p.value = PVAL, method = METHOD,
                 data.name = DNAME, parameter = PARAMETER), 
            class = "bdstest")
}

print.bdstest <- function (x, digits = 4)
{
  cat("\n\t", x$method, "\n\n")
  cat("data: ", x$data.name, "\n\n")
  if (!is.null(x$parameter))
  {
    cat("Embedding dimension = ", format(round(x$parameter$m, digits)), sep = " ", "\n\n")
    cat("Epsilon for close points = ", format(round(x$parameter$eps, digits)), sep = " ", "\n\n")
  }
  if (!is.null(x$statistic))
  {
    colnames(x$statistic) <- round (as.numeric(colnames(x$statistic)), digits)
    colnames(x$statistic) <- paste("[",colnames(x$statistic),"]")
    rownames(x$statistic) <- round (as.numeric(rownames(x$statistic)), digits)
    rownames(x$statistic) <- paste("[",rownames(x$statistic),"]")
    cat("Standard Normal = \n")
    print (round(x$statistic, digits))
    cat("\n")
  }
  if (!is.null(x$p.value))
  {
    colnames(x$p.value) <- round (as.numeric(colnames(x$p.value)), digits)
    colnames(x$p.value) <- paste("[",colnames(x$p.value),"]")
    rownames(x$p.value) <- round (as.numeric(rownames(x$p.value)), digits)
    rownames(x$p.value) <- paste("[",rownames(x$p.value),"]")
    cat("p-value = \n")
    print (round(x$p.value, digits))
    cat("\n")
  }
  cat("\n")
  invisible(x)
}

adf.test <- function (x, k = trunc((length(x)-1)^(1/3)))
{
  if (!is.vector(x) & !is.univariate.ts(x))
    stop ("x is not a vector or univariate time series")
  if (k < 0) stop ("k negative")
  DNAME <- deparse(substitute(x))
  k <- k+1
  y <- diff (x)
  n <- length(y)
  z <- embed (y, k)
  yt <- z[,1]
  xt1 <- x[k:n]
  tt <- k:n
  if (k > 1)
  {
    yt1 <- z[,2:k] 
    res <- lm (yt~xt1+1+tt+yt1)
  }
  else
    res <- lm (yt~xt1+1+tt)
  res.sum <- summary (res)
  STAT <- res.sum$coefficients[2,1]/res.sum$coefficients[2,2]
  table <- cbind(c(4.38,4.15,4.04,3.99,3.98,3.96),
                 c(3.95,3.80,3.73,3.69,3.68,3.66),
                 c(3.60,3.50,3.45,3.43,3.42,3.41),
                 c(3.24,3.18,3.15,3.13,3.13,3.12),
                 c(1.14,1.19,1.22,1.23,1.24,1.25),
                 c(0.80,0.87,0.90,0.92,0.93,0.94),
                 c(0.50,0.58,0.62,0.64,0.65,0.66),
                 c(0.15,0.24,0.28,0.31,0.32,0.33))
  table <- -table
  tablen <- dim(table)[2]
  tableT <- c(25,50,100,250,500,100000)
  tablep <- c(0.01,0.025,0.05,0.10,0.90,0.95,0.975,0.99)
  tableipl <- numeric(tablen)
  for (i in (1:tablen))
    tableipl[i] <- approx (tableT,table[,i],n,rule=2)$y
  PVAL <- approx (tableipl,tablep,STAT,rule=2)$y
  PARAMETER <- k-1
  METHOD <- "Augmented Dickey-Fuller Test"
  names(STAT) <- "Dickey-Fuller"
  names(PARAMETER) <- "Lag order"
  structure(list(statistic = STAT, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name = DNAME), 
            class = "htest")
}
