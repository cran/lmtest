bptest <- function(formula, varformula = NULL, studentize = TRUE,
  data = list(), weights = NULL)
{
  dname <- paste(deparse(substitute(formula)))

  if(!inherits(formula, "formula")) {
    X <- if(is.matrix(formula$x))
	   formula$x
	 else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
	   formula$y
	 else model.response(model.frame(formula))
    Z <- if(is.null(varformula)) X
           else model.matrix(varformula, data = data)
    wts <- weights(formula)
  } else {
    mf <- if(is.null(weights)) {
      model.frame(formula, data = data)
    } else {
      model.frame(formula, weights = weights, data = data)
    }
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
    Z <- if(is.null(varformula)) X
           else model.matrix(varformula, data = data)
    wts <- model.weights(mf)
  }
  if(is.null(wts)) wts <- rep.int(1, NROW(X))

  ## only use complete cases that are in both models
  if(!(all(c(row.names(X) %in% row.names(Z), row.names(Z) %in% row.names(X))))) {
    allnames <- row.names(X)[row.names(X) %in% row.names(Z)]
    X <- X[allnames,]
    Z <- Z[allnames,]
    y <- y[allnames]
    wts <- wts[row.names(X) %in% row.names(Z)]
  }

  ## need at least one intercept plus one regressor
  if(ncol(Z) < 2L) stop("the auxiliary variance regression requires at least an intercept and a regressor")
   
  k <- ncol(X)
  n <- sum(wts > 0) ## nrow(X)

  resi <- lm.wfit(X, y, wts)$residuals
  sigma2 <- sum(wts * resi^2)/n

  if(studentize) {
    w <- resi^2 - sigma2
    aux <- lm.wfit(Z, w, wts)    
    bp <- n * sum(wts * aux$fitted.values^2)/sum(w^2)
    method <- "studentized Breusch-Pagan test"
  } else {
    f <- resi^2/sigma2 -1
    aux <- lm.wfit(Z, f, wts)
    bp <- 0.5 * sum(wts * aux$fitted.values^2)
    method <- "Breusch-Pagan test"
  }

  names(bp) <- "BP"
  df <- c("df" = aux$rank - 1)
  RVAL <- list(statistic = bp,
      parameter = df,
      method = method,
      p.value= pchisq(bp, df, lower.tail = FALSE),
      data.name = dname)

  class(RVAL) <- "htest"
  return(RVAL)
}

