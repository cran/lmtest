reset <- resettest <- function(formula, power = 2:3,
  type = c("fitted", "regressor", "princomp"), data = list(),
  vcov = NULL, ...)
{
  dname <- paste(deparse(substitute(formula)))

  if(!inherits(formula, "formula")) {
    mf <- model.frame(formula)
    X <- if(is.matrix(formula$x))
           formula$x
         else model.matrix(terms(formula), mf)
    y <- if(is.vector(formula$y))
           formula$y
         else model.response(mf)
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }
  y <- scale(y, center = attr(terms(mf), "intercept") > 0L)
  k <- ncol(X)
  n <- nrow(X)
  type <- match.arg(type)
  switch(type,

  "fitted" = {
    y.hat <- lm.fit(X,y)$fitted.values
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

  if(is.null(vcov)) {
    XZ <- cbind(X, Z)
    q <- ncol(Z)
    res1 <- lm.fit(X,y)$residuals
    res2 <- lm.fit(XZ,y)$residuals
    res1 <- sum(res1^2)
    res2 <- sum(res2^2)
    df1 <- q
    df2 <- n-(k+q)
    reset <- (df2/df1) * ((res1 - res2) / res2)
    df <- c(df1, df2)
    pval <- as.vector(pf(reset, df1, df2, lower.tail = FALSE))
  } else {
    reset <- waldtest(lm(y ~ 0 + X + Z), "Z", vcov = vcov, ...)
    df <- abs(as.numeric(reset[2L, 2L:1L]))
    pval <- as.numeric(reset[2L, 4L])
    reset <- as.numeric(reset[2L, 3L])
  }
  names(reset) <- "RESET"
  names(df) <- c("df1", "df2")
  
  RVAL <- list(statistic = reset,
	parameter = df,
	method = "RESET test",
	p.value = pval,
	data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
}
