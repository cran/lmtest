coeftest <- function(x, vcov = NULL, df = NULL, ...)
{
  UseMethod("coeftest")
}

coeftest.default <- function(x, vcov = NULL, df = NULL, ...)
{
  est <- as.vector(coef(x))
  if(is.null(vcov)) se <- vcov(x) else {
      if(is.function(vcov)) se <- vcov(x)
        else se <- vcov
  }
  se <- sqrt(diag(se))
  tval <- est/se

  if(is.null(df)) df <- x$df.residual
  if(is.null(df)) df <- 0

  if(is.finite(df) && df > 0) {
    pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    mthd <- "t"
  } else {
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    mthd <- "z"
  }
  rval <- cbind(est, se, tval, pval)
  colnames(rval) <- cnames
  class(rval) <- "coeftest"
  attr(rval, "method") <- paste(mthd, "test of coefficients of",
    dQuote(class(x)[1]), "object", sQuote(deparse(substitute(x))))
  return(rval)
} 

coeftest.breakpointsfull <- function(x, vcov = NULL, df = NULL, ...)
{
  stopifnot(require(strucchange))
  
  est <- coef(x, ...)
  rnames <- as.vector(t(outer(rownames(est), colnames(est), paste)))
  est <- as.vector(t(est))
  
  se <- vcov.breakpointsfull(x, vcov = vcov, ...)

  se <- as.vector(sapply(seq(along = se), function(x) sqrt(diag(se[[x]]))))
  tval <- est/se

  if(is.finite(df) && df > 0) {
    pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    mthd <- "t"
  } else {
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)
    cnames <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    mthd <- "z"
  }
  rval <- cbind(est, se, tval, pval)
  colnames(rval) <- cnames
  rownames(rval) <- rnames
  class(rval) <- "coeftest"
  attr(rval, "method") <- paste(mthd, "test of coefficients of",
    dQuote(class(x)[1]), "object", sQuote(deparse(substitute(x))))
  return(rval)
} 

print.coeftest <- function(x, ...)
{
  mthd <- attr(x, "method")
  if(is.null(mthd)) mthd <- "Test of coefficients"
  cat(paste("\n", mthd,":\n\n", sep = ""))
  printCoefmat(x, ...)
  cat("\n")
  invisible(x)
}
