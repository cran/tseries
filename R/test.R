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
# Mostly time series tests
#


runs.test <- function (x)
{
  if (NCOL(x) > 1) stop ("x is not a vector or univariate time series")
  if (any(is.na(x))) stop ("NAs in x")
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
  METHOD <- "Runs Test"
  PVAL <- 2 * pnorm (-abs(STATISTIC))
  names(STATISTIC) <- "Standard Normal"
  structure(list(statistic = STATISTIC,
		 p.value = PVAL,
		 method = METHOD,
		 data.name = DNAME),
	    class = "htest")
}

bds.test <- function (x, m = 2, eps = seq(0.5*sd(x),2*sd(x),length=4), trace = FALSE)
{
  if (NCOL(x) > 1) stop ("x is not a vector or univariate time series")
  if (any(is.na(x))) stop ("NAs in x")
  if (m < 2) stop ("m is less than 2")
  if (length(eps) == 0) stop ("invalid eps")
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
              as.vector(cc), cstan=as.vector(cstan), as.double(eps[i]), as.integer(trace),
              PACKAGE="tseries")
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

print.bdstest <- function (object, digits = 4)
{
  if (!inherits(object, "bdstest")) stop ("method is only for bdstest objects")
  cat("\n\t", object$method, "\n\n")
  cat("data: ", object$data.name, "\n\n")
  if (!is.null(object$parameter))
  {
    cat("Embedding dimension = ", format(round(object$parameter$m, digits)), sep = " ", "\n\n")
    cat("Epsilon for close points = ", format(round(object$parameter$eps, digits)),
        sep = " ", "\n\n")
  }
  if (!is.null(object$statistic))
  {
    colnames(object$statistic) <- round (as.numeric(colnames(object$statistic)), digits)
    colnames(object$statistic) <- paste("[",colnames(object$statistic),"]")
    rownames(object$statistic) <- round (as.numeric(rownames(object$statistic)), digits)
    rownames(object$statistic) <- paste("[",rownames(object$statistic),"]")
    cat("Standard Normal = \n")
    print (round(object$statistic, digits))
    cat("\n")
  }
  if (!is.null(object$p.value))
  {
    colnames(object$p.value) <- round (as.numeric(colnames(object$p.value)), digits)
    colnames(object$p.value) <- paste("[",colnames(object$p.value),"]")
    rownames(object$p.value) <- round (as.numeric(rownames(object$p.value)), digits)
    rownames(object$p.value) <- paste("[",rownames(object$p.value),"]")
    cat("p-value = \n")
    print (round(object$p.value, digits))
    cat("\n")
  }
  cat("\n")
  invisible(object)
}

adf.test <- function (x, k = trunc((length(x)-1)^(1/3)))
{
  if (NCOL(x) > 1) stop ("x is not a vector or univariate time series")
  if (any(is.na(x))) stop ("NAs in x")
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

white.test <- function (object, ...) { UseMethod("white.test") }

white.test.default <- function (x, y, qstar = 2, q = 10, range = 4,
                                type = c("chisq","F"), scale = TRUE)
{
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (any(is.na(x))) stop ("NAs in x")
  if (any(is.na(y))) stop ("NAs in y")
  nin <- dim(x)[2]
  t <- dim(x)[1]
  if (dim(x)[1] != dim(y)[1]) 
    stop("number of rows of x and y must match")
  if (dim(x)[1] <= 0) 
    stop("no observations in x and y")
  if (dim(y)[2] > 1)
    stop ("handles only univariate outputs")
  if (!require (mva, quietly=TRUE))
    stop ("Package mva is needed. Stopping")
  type <- match.arg (type)
  if (scale)
  {
    x <- scale(x)
    y <- scale(y)
  }
  xnam <- paste ("x[,", 1:nin, "]", sep="")
  fmla <- as.formula (paste ("y~",paste(xnam,collapse= "+")))
  rr <- lm (fmla)
  u <- residuals (rr)
  ssr0 <- sum (u^2)
  max <- range/2
  gamma <- matrix(runif((nin+1)*q,-max,max),nin+1,q)
  phantom <- (1+exp(-(cbind(rep(1,t),x)%*%gamma)))^(-1)
  phantomstar <- as.matrix(prcomp(phantom,scale=TRUE)$x[,2:(qstar+1)])
  xnam2 <- paste ("phantomstar[,", 1:qstar, "]", sep="")
  xnam2 <- paste(xnam2,collapse="+")
  fmla <- as.formula (paste ("u~",paste(paste(xnam,collapse= "+"),xnam2,sep="+")))
  rr <- lm (fmla)
  v <- residuals(rr)
  ssr <- sum(v^2)
  if (type == "chisq")
  {
    STAT <- t*log(ssr0/ssr)
    PVAL <- 1-pchisq(STAT,qstar)
    PARAMETER <- qstar
    names(STAT) <- "X-squared"
    names(PARAMETER) <- "df"
  }
  else if (type == "F")
  {
    STAT <- ((ssr0-ssr)/qstar)/(ssr/(t-qstar-nin))
    PVAL <- 1-pf(STAT,qstar,t-qstar-nin)
    PARAMETER <- c(qstar,t-qstar-nin)
    names(STAT) <- "F"
    names(PARAMETER) <- c("df1","df2")
  }
  else
    stop ("invalid type")
  ARG <- c(qstar,q,range,scale)
  names(ARG) <- c("qstar","q","range","scale")
  METHOD <- "White Neural Network Test"
  structure(list(statistic = STAT, parameter = PARAMETER, p.value = PVAL, 
                 method = METHOD, data.name = DNAME, arguments = ARG), class = "htest")
}

white.test.ts <- function (x, lag = 1, qstar = 2, q = 10, range = 4,
                           type = c("chisq","F"), scale = TRUE)
{
  if (!is.ts(x)) stop ("method is only for time series")
  if (NCOL(x) > 1) stop ("x is not a vector or univariate time series")
  if (any(is.na(x))) stop ("NAs in x")
  if (lag < 1) 
    stop("minimum lag is 1")
  if (!require (mva, quietly=TRUE))
    stop ("Package mva is needed. Stopping")
  type <- match.arg (type)
  DNAME <- deparse(substitute(x))
  t <- length(x)
  if (scale) x <- scale(x)
  y <- embed (x, lag+1)
  xnam <- paste ("y[,", 2:(lag+1), "]", sep="")
  fmla <- as.formula (paste ("y[,1]~",paste(xnam,collapse= "+")))
  rr <- lm (fmla)
  u <- residuals (rr)
  ssr0 <- sum (u^2)
  max <- range/2
  gamma <- matrix(runif((lag+1)*q,-max,max),lag+1,q)
  phantom <- (1+exp(-(cbind(rep(1,t-lag),y[,2:(lag+1)])%*%gamma)))^(-1)
  phantomstar <- as.matrix(prcomp(phantom,scale=TRUE)$x[,2:(qstar+1)])
  xnam2 <- paste ("phantomstar[,", 1:qstar, "]", sep="")
  xnam2 <- paste(xnam2,collapse="+")
  fmla <- as.formula (paste ("u~",paste(paste(xnam,collapse= "+"),xnam2,sep="+")))
  rr <- lm (fmla)
  v <- residuals(rr)
  ssr <- sum(v^2)
  if (type == "chisq")
  {
    STAT <- t*log(ssr0/ssr)
    PVAL <- 1-pchisq(STAT,qstar)
    PARAMETER <- qstar
    names(STAT) <- "X-squared"
    names(PARAMETER) <- "df"
  }
  else if (type == "F")
  {
    STAT <- ((ssr0-ssr)/qstar)/(ssr/(t-lag-qstar))
    PVAL <- 1-pf(STAT,qstar,t-lag-qstar)
    PARAMETER <- c(qstar,t-lag-qstar)
    names(STAT) <- "F"
    names(PARAMETER) <- c("df1","df2")
  }
  else
    stop ("invalid type")
  ARG <- c(lag,qstar,q,range,scale)
  names(ARG) <- c("lag","qstar","q","range","scale")
  METHOD <- "White Neural Network Test"
  structure(list(statistic = STAT, parameter = PARAMETER, p.value = PVAL, 
                 method = METHOD, data.name = DNAME, arguments = ARG), class = "htest")
}

terasvirta.test <- function (object, ...) { UseMethod("terasvirta.test") }

terasvirta.test.default <- function (x, y, type = c("chisq","F"), scale = TRUE)
{
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  x <- as.matrix(x)
  y <- as.matrix(y)
  if (any(is.na(x))) stop ("NAs in x")
  if (any(is.na(y))) stop ("NAs in y")
  nin <- dim(x)[2]
  if (nin < 1)
    stop ("invalid x")
  t <- dim(x)[1]
  if (dim(x)[1] != dim(y)[1]) 
    stop("number of rows of x and y must match")
  if (dim(x)[1] <= 0) 
    stop("no observations in x and y")
  if (dim(y)[2] > 1)
    stop ("handles only univariate outputs")
  type <- match.arg (type)
  if (scale)
  {
    x <- scale(x)
    y <- scale(y)
  }
  xnam <- paste ("x[,", 1:nin, "]", sep="")
  fmla <- as.formula (paste ("y~",paste(xnam,collapse= "+")))
  rr <- lm (fmla)
  u <- residuals (rr)
  ssr0 <- sum (u^2)
  xnam2 <- NULL
  m <- 0
  for (i in (1:nin))
  {
    for (j in (i:nin))
    {
      xnam2 <- c(xnam2,paste("I(x[,",i,"]*x[,",j,"])",sep=""))
      m <- m+1
    }
  }
  xnam2 <- paste(xnam2,collapse="+")
  xnam3 <- NULL
  for (i in (1:nin))
  {
    for (j in (i:nin))
    {
      for (k in (j:nin))
      {
        xnam3 <- c(xnam3,paste("I(x[,",i,"]*x[,",j,"]*x[,",k,"])",sep=""))
        m <- m+1
      }
    }
  }
  xnam3 <- paste(xnam3,collapse="+")
  fmla <- as.formula (paste ("u~",paste(paste(xnam,collapse= "+"),xnam2,xnam3,sep="+")))
  rr <- lm (fmla)
  v <- residuals(rr)
  ssr <- sum(v^2)
  if (type == "chisq")
  {
    STAT <- t*log(ssr0/ssr)
    PVAL <- 1-pchisq(STAT,m)
    PARAMETER <- m
    names(STAT) <- "X-squared"
    names(PARAMETER) <- "df"
  }
  else if (type == "F")
  {
    STAT <- ((ssr0-ssr)/m)/(ssr/(t-nin-m))
    PVAL <- 1-pf(STAT,m,t-nin-m)
    PARAMETER <- c(m,t-nin-m)
    names(STAT) <- "F"
    names(PARAMETER) <- c("df1","df2")
  }
  else
    stop ("invalid type")
  METHOD <- "Teraesvirta Neural Network Test"
  ARG <- scale
  names(ARG) <- "scale"
  structure(list(statistic = STAT, parameter = PARAMETER, p.value = PVAL, 
                 method = METHOD, data.name = DNAME, arguments = ARG), class = "htest")
}

terasvirta.test.ts <- function (x, lag = 1, type = c("chisq","F"), scale = TRUE)
{
  if (!is.ts(x)) stop ("method is only for time series")
  if (NCOL(x) > 1) stop ("x is not a vector or univariate time series")
  if (any(is.na(x))) stop ("NAs in x")
  if (lag < 1) 
    stop("minimum lag is 1")
  type <- match.arg (type)
  DNAME <- deparse(substitute(x))
  t <- length(x)
  if (scale) x <- scale(x)
  y <- embed (x, lag+1)
  xnam <- paste ("y[,", 2:(lag+1), "]", sep="")
  fmla <- as.formula (paste ("y[,1]~",paste(xnam,collapse= "+")))
  rr <- lm (fmla)
  u <- residuals (rr)
  ssr0 <- sum (u^2)
  xnam2 <- NULL
  m <- 0
  for (i in (1:lag))
  {
    for (j in (i:lag))
    {
      xnam2 <- c(xnam2,paste("I(y[,",i+1,"]*y[,",j+1,"])",sep=""))
      m <- m+1
    }
  }
  xnam2 <- paste(xnam2,collapse="+")
  xnam3 <- NULL
  for (i in (1:lag))
  {
    for (j in (i:lag))
    {
      for (k in (j:lag))
      {
        xnam3 <- c(xnam3,paste("I(y[,",i+1,"]*y[,",j+1,"]*y[,",k+1,"])",sep=""))
        m <- m+1
      }
    }
  }
  xnam3 <- paste(xnam3,collapse="+")
  fmla <- as.formula (paste ("u~",paste(paste(xnam,collapse= "+"),xnam2,xnam3,sep="+")))
  rr <- lm (fmla)
  v <- residuals(rr)
  ssr <- sum(v^2)
  if (type == "chisq")
  {
    STAT <- t*log(ssr0/ssr)
    PVAL <- 1-pchisq(STAT,m)
    PARAMETER <- m
    names(STAT) <- "X-squared"
    names(PARAMETER) <- "df"
  }
  else if (type == "F")
  {
    STAT <- ((ssr0-ssr)/m)/(ssr/(t-lag-m))
    PVAL <- 1-pf(STAT,m,t-lag-m)
    PARAMETER <- c(m,t-lag-m)
    names(STAT) <- "F"
    names(PARAMETER) <- c("df1","df2")
  }
  else
    stop ("invalid type")
  METHOD <- "Teraesvirta Neural Network Test"
  ARG <- c(lag,scale)
  names(ARG) <- c("lag","scale")
  structure(list(statistic = STAT, parameter = PARAMETER, p.value = PVAL, 
                 method = METHOD, data.name = DNAME, arguments = ARG), class = "htest")
}

jarque.bera.test <- function (x)
{
  if (NCOL(x) > 1) stop ("x is not a vector or univariate time series")
  if (any(is.na(x))) stop ("NAs in x")
  DNAME <- deparse (substitute(x))
  n <- length (x)
  m1 <- sum(x)/n
  m2 <- sum((x-m1)^2)/n
  m3 <- sum((x-m1)^3)/n
  m4 <- sum((x-m1)^4)/n
  b1 <- (m3/m2^(3/2))^2
  b2 <- (m4/m2^2)
  STATISTIC <- n*b1/6+n*(b2-3)^2/24
  names(STATISTIC) <- "X-squared"
  PARAMETER <- 2
  names(PARAMETER) <- "df"
  PVAL <- 1-pchisq(STATISTIC,df = 2)
  METHOD <- "Jarque Bera Test"
  structure(list(statistic = STATISTIC,
                 parameter = PARAMETER,
		 p.value = PVAL,
		 method = METHOD,
		 data.name = DNAME),
	    class = "htest")
}

