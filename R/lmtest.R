dwtest <- function(formula, alternative = c("greater", "two.sided", "less"),
  iterations = 15, exact = NULL, tol = 1e-10, data = list())
{
  dname <- paste(deparse(substitute(formula)))
  alternative <- match.arg(alternative)

  if(!"formula" %in% class(formula)) {
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  
   
  n <- nrow(X)
  if(is.null(exact)) exact <- (n < 100)
  k <- ncol(X)
  res <- lm.fit(X,y)$residuals
  dw <- sum(diff(res)^2)/sum(res^2)
  A <- diag(c(1,rep(2, n-2), 1))
  A[abs(row(A)-col(A))==1] <- -1
  Q1 <- chol2inv(qr.R(qr(X)))
  if(exact)
  {
    MA <- diag(rep(1,n)) - X %*% Q1 %*% t(X)
    MA <- MA %*% A
    ev <- eigen(MA)$values[1:(n-k)]
    if(any(Im(ev)>tol)) warning("imaginary parts of eigenvalues discarded")
    ev <- Re(ev)
    ev <- ev[ev>tol]
    pdw <- function(dw) .Fortran("pan", as.double(c(dw,ev)), as.integer(length(ev)),
             as.double(0), as.integer(iterations), x=double(1), PACKAGE = "lmtest")$x
    pval <- switch(alternative,
      "two.sided" = (2*min(pdw(dw), 1-pdw(dw))),
      "less" = (1 - pdw(dw)),
      "greater" = pdw(dw))

    if((pval > 1) | (pval < 0))
    {
      warning("exact p value cannot be computed (not in [0,1]), approximate p value will be used")
      exact <- FALSE
    }
  }
  if(!exact)
  {
    XAXQ <- t(X) %*% A %*% X %*% Q1
    P <- 2*(n-1) - sum(diag(XAXQ))
    Q <- 2*(3*n - 4) - 2* sum(diag(t(X) %*% A %*% A %*% X %*% Q1)) + sum(
         diag(XAXQ %*% XAXQ))
    dmean <- P/(n-k)
    dvar <- 2/((n-k)*(n-k+2)) * (Q - P*dmean)
    pval <- switch(alternative,
      "two.sided" = (2*pnorm(abs(dw-dmean), sd=sqrt(dvar), lower.tail = FALSE)),
      "less" = pnorm(dw, mean = dmean, sd = sqrt(dvar), lower.tail = FALSE),
      "greater" = pnorm(dw, mean = dmean, sd = sqrt(dvar)))
  }

  alternative <- switch(alternative,
    "two.sided" = "true autocorelation is not 0",
    "less" = "true autocorrelation is less than 0",
    "greater" = "true autocorrelation is greater than 0")

  names(dw) <- "DW"
  RVAL <- list(statistic = dw, method = "Durbin-Watson test",
    alternative = alternative, p.value= pval, data.name=dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

bptest <- function(formula, varformula = NULL, studentize = TRUE,
  data = list())
{
  dname <- paste(deparse(substitute(formula)))

  if(!"formula" %in% class(formula)) {
    X <- if(is.matrix(formula$x))
	   formula$x
	 else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
	   formula$y
	 else model.response(model.frame(formula))
    Z <- if(is.null(varformula)) X
           else model.matrix(varformula, data = data)
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    Z <- if(is.null(varformula)) X
           else model.matrix(varformula, data = data)
  }  
   
  k <- ncol(X)
  n <- nrow(X)

  resi <- lm.fit(X,y)$residuals
  sigma2 <- sum(resi^2)/n

  if(studentize)
  {
    w <- resi^2 - sigma2
    fv <- lm.fit(Z,w)$fitted
    bp <- n * sum(fv^2)/sum(w^2)
    method <- "studentized Breusch-Pagan test"
  }
  else
  {
    f <- resi^2/sigma2 -1
    fv <- lm.fit(Z,f)$fitted
    bp <- 0.5 * sum(fv^2)
    method <- "Breusch-Pagan test"
  }

  names(bp) <- "BP"
  df <- ncol(Z)-1
  names(df) <- "df";
  RVAL <- list(statistic = bp,
      parameter = df,
      method = method,
      p.value= 1-pchisq(bp,df),
      data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

gqtest <- function(formula, point = 0.5, order.by = NULL, data = list())
{
  dname <- paste(deparse(substitute(formula)))

  if(!"formula" %in% class(formula)) {
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  

  k <- ncol(X)
  n <- nrow(X)
  if(point < 1) point <- floor(point*n)
  if (point > n - k | point < k) stop("inadmissable breakpoint")

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
  }

  rss1 <- sum(lm.fit(as.matrix(X[1:point,]),y[1:point])$residuals^2)
  rss2 <- sum(lm.fit(as.matrix(X[(point+1):n,]),y[(point+1):n])$residuals^2)

  gq <- (rss2/(n-point-k))/(rss1/(point-k))
  df <- c(n-point-k, point-k)
  names(df) <- c("df1", "df2")
  PVAL <- 1-pf(gq, df[1], df[2])
  method <- "Goldfeld-Quandt test"
  names(gq) <- "GQ"
  RVAL <- list(statistic = gq,
      parameter = df,
      method = method,
      p.value= PVAL,
      data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

hmctest <- function(formula, point = 0.5, order.by = NULL, simulate.p = TRUE,
  nsim = 1000, plot = FALSE, data = list())
{
  dname <- paste(deparse(substitute(formula)))
  
  if(!"formula" %in% class(formula)) {
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  

  k <- ncol(X)
  n <- nrow(X)
  if(point < 1) point <- floor(point*n)
  if (point > n - k | point < k) stop("inadmissable breakpoint")

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
  }


  resi <- lm.fit(X,y)$residuals
  hmc <- sum(resi[1:point]^2)/sum(resi^2)

  if(plot)
  {
    stats <- c(0,cumsum(resi^2))/sum(resi^2)
    stats <- ts(stats, start=0, freq=n)
    plot(stats, xlab="fraction", ylab="Harrison-McCabe statistics", xaxs="i",
      yaxs="i")
    abline(0,1)
  }

  names(hmc) <- "HMC"
  if (simulate.p)
  {
    stat <- rep(0, nsim)
    for (i in 1:nsim) {
      x <- rnorm(n)
      x <- (x - mean(x))/sqrt(var(x))
      stat[i] <- sum(x[1:point]^2)/sum(x^2)
    }
    PVAL <- mean(stat <= hmc)
  }
  else
    PVAL <- NA

  RVAL <- list(statistic = hmc,
      method = "Harrison-McCabe test",
      p.value= PVAL,
      data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

harvtest <- function(formula, order.by=NULL, data=list())
{
  dname <- paste(deparse(substitute(formula)))

  if(!"formula" %in% class(formula)) {
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  

  k <- ncol(X)
  n <- nrow(X)

  rec.res <- function(X, y)
  {
      n <- nrow(X)
      q <- ncol(X)
      w <- rep(0,(n-q))
      Xr1 <- X[1:q,,drop = FALSE]
      xr <- as.vector(X[q+1,])
      X1 <- chol2inv(qr.R(qr(Xr1)))
      fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
      betar <- X1 %*%t(Xr1)%*% y[1:q]
      w[1] <- (y[q+1] - t(xr) %*% betar)/sqrt(fr)

      for(r in ((q+2):n))
      {
          X1 <- X1 - (X1 %*% outer(xr, xr) %*% X1)/fr
  	  betar <- betar + X1 %*% xr * w[r-q-1]*sqrt(fr)
  	  xr <- as.vector(X[r,])
          fr <- as.vector((1 + (t(xr) %*% X1 %*% xr)))
          w[r-q] <- (y[r] - t(xr) %*% betar)/sqrt(fr)
      }
      return(w)
  }

  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
  }

  resr <- rec.res(X,y)
  sigma <- sqrt(var(resr)*(length(resr)-1)/(n-k-1))
  resr <- resr / sigma
  harv <- abs(sum(resr)/sqrt(n-k))/sqrt(var(resr))
  names(harv) <- "HC"
  df <- n-k-1
  names(df) <- "df"
  RVAL <- list(statistic = harv,
      parameter = df,
      method = "Harvey-Collier test",
      p.value= 2 * (1-pt(harv, n-k-1)),
      data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

raintest <- function(formula, fraction=0.5, order.by=NULL, center=NULL,
 data=list())
{
  dname <- paste(deparse(substitute(formula)))

  if(!"formula" %in% class(formula)) {
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  

  k <- ncol(X)
  n <- nrow(X)

  if(is.null(order.by))
  {
    if(is.null(center)) center <- 0.5
    if(center > 1) center <- center/n
    from <- ceiling(quantile(1:n, probs=(center-fraction/2)))
    to <- from + floor(fraction*n) - 1
  }
  else
  if(order.by == "mahalanobis")
  {
    if(is.null(center)) center <- apply(X,2,mean)
    o <- order(mahalanobis(X,center, chol2inv(qr.R(qr(X))), inverted = TRUE))
    X <- as.matrix(X[o,])
    y <- y[o]
    from <- 1
    to <- floor(fraction*n)
  }
  else
  {
    if(inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[,ncol(z)])
    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
    if(is.null(center)) center <- 0.5
    if(center > 1) center <- center/n
    from <- ceiling(quantile(1:n, probs=(center-fraction/2)))
    to <- from + floor(fraction*n) - 1
  }

  subX <- as.matrix(X[from:to,])
  suby <- y[from:to]
  n1 <- nrow(subX)
  if(n1 < k) stop("not enough observations in subset")
  resi <- lm.fit(X,y)$residuals
  subresi <- lm.fit(subX, suby)$residuals
  sresi <- sum(resi^2)
  sresi1 <- sum(subresi^2)
  rain <- ((sresi - sresi1)/(n-n1))/(sresi1/(n1-k))
  names(rain) <- "Rain"
  df <- c((n-n1),(n1-k))
  names(df) <- c("df1","df2")
  RVAL <- list(statistic = rain,
      parameter = df,
      method = "Rainbow test",
      p.value= as.vector(1-pf(rain, df[1], df[2])),
      data.name=dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

reset <- function(formula, power=2:3, type=c("fitted", "regressor",
  "princomp"), data=list())
{
  dname <- paste(deparse(substitute(formula)))

  if(!"formula" %in% class(formula)) {
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  

  k <- ncol(X)
  n <- nrow(X)
  type <- match.arg(type)
  switch(type,

  "fitted" = {
    y.hat <- lm.fit(X,y)$fitted
    Z <- matrix(t(sapply(y.hat, "^", power)), nrow=n)
  },

  "regressor" = {
    Z <- as.matrix(mf[,which(!sapply(mf,is.factor))[-1]])
    Z <- matrix(as.vector(t(sapply(as.vector(Z), "^", power))), nrow=n)
  },

  "princomp" = {
    Z <- as.matrix(mf[,which(!sapply(mf,is.factor))[-1]])
    pc1 <- as.matrix(eigen(cov(Z))$vectors)[,1]
    pc1 <- as.vector(Z %*% pc1)
    Z <- matrix(t(sapply(pc1, "^", power)), nrow=n)
  })

  XZ <- cbind(X, Z)
  q <- ncol(Z)
  res1 <- lm.fit(X,y)$residuals
  res2 <- lm.fit(XZ,y)$residuals
  res1 <- sum(res1^2)
  res2 <- sum(res2^2)
  df1 <- q
  df2 <- n-(k+q)
  reset <- (df2/df1) * ((res1 - res2) / res2)
  names(reset) <- "RESET"
  df <- c(df1, df2)
  names(df) <- c("df1","df2")
  RVAL <- list(statistic = reset,
	parameter = df,
	method = "RESET test",
	p.value= as.vector(1-pf(reset, df1, df2)),
	data.name=dname)
  class(RVAL) <- "htest"
  return(RVAL)
}

bgtest <- function(formula, order = 1, type = c("Chisq", "F"), data = list())
{
  dname <- paste(deparse(substitute(formula)))

  if(!"formula" %in% class(formula)) {
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  

  n <- nrow(X)
  k <- ncol(X)
  order <- 1:order
  m <- length(order)
  resi <- lm.fit(X,y)$residuals

  Z <- sapply(order, function(x) c(rep(0, x), resi[1:(n-x)]))
  auxfit <- lm.fit(cbind(X,Z), resi)

  switch(match.arg(type),

  "Chisq" = {
    bg <- n * sum(auxfit$fitted^2)/sum(resi^2)
    p.val <- pchisq(bg, m, lower.tail = FALSE)
    df <- m
    names(df) <- "df"
  },

  "F" = {
    uresi <- auxfit$residuals
    rresi <- lm.fit(X,resi)$residuals
    bg <- ((sum(rresi^2) - sum(uresi^2))/m) / (sum(uresi^2) / (n-k-m))
    p.val <- pf(bg, df1 = m, df2 = n-k, lower.tail = FALSE)
    df <- c(m, n-k)
    names(df) <- c("df1", "df2")
  })

  names(bg) <- "LM test"
  RVAL <- list(statistic = bg, parameter = df,
               method = paste("Breusch-Godfrey test for serial correlation of order", max(order)),
               p.value = p.val,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)
}
