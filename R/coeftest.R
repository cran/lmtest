coeftest <- function(x, vcov. = NULL, df = NULL, ...)
{
  UseMethod("coeftest")
}

coeftest.default <- function(x, vcov. = NULL, df = NULL, ..., save = FALSE)
{
  ## use S4 methods if loaded
  coef0 <- if("stats4" %in% loadedNamespaces()) stats4::coef else coef
  vcov0 <- if("stats4" %in% loadedNamespaces()) stats4::vcov else vcov
  nobs0 <- if("stats4" %in% loadedNamespaces()) stats4::nobs else nobs
  logl0 <- if("stats4" %in% loadedNamespaces()) stats4::logLik else logLik

  ## extract coefficients and standard errors
  est <- coef0(x)
  if(is.null(vcov.)) se <- vcov0(x) else {
      if(is.function(vcov.)) se <- vcov.(x, ...)
        else se <- vcov.
  }
  se <- sqrt(diag(se))

  ## match using names and compute t/z statistics
  if(!is.null(names(est)) && !is.null(names(se))) {
    if(length(unique(names(est))) == length(names(est)) && length(unique(names(se))) == length(names(se))) {
      anames <- names(est)[names(est) %in% names(se)]
      est <- est[anames]
      se <- se[anames]
    }
  }  
  tval <- as.vector(est)/se

  ## apply central limit theorem
  if(is.null(df)) {
    df <- try(df.residual(x), silent = TRUE)
    if(inherits(df, "try-error")) df <- NULL
  }
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
  attr(rval, "method") <- paste(mthd, "test of coefficients")
  attr(rval, "df") <- df
  
  ## supplementary information for model summary
  n <- try(nobs0(x), silent = TRUE)
  attr(rval, "nobs") <- if(inherits(n, "try-error")) NULL else n
  ll <- try(logl0(x), silent = TRUE)
  attr(rval, "logLik") <- if(inherits(ll, "try-error")) NULL else ll
  if(save) attr(rval, "object") <- x
  
  return(rval)
} 

coeftest.glm <- function(x, vcov. = NULL, df = Inf, ...)
  coeftest.default(x, vcov. = vcov., df = df, ...)  

coeftest.mlm <- function(x, vcov. = NULL, df = NULL, ...)
{
  ## obtain vcov
  v <- if(is.null(vcov.)) vcov(x) else if(is.function(vcov.)) vcov.(x) else vcov.

  ## nasty hack: replace coefficients so that their names match the vcov() method
  x$coefficients <- structure(as.vector(x$coefficients), .Names = colnames(vcov(x)))

  ## call default method
  coeftest.default(x, vcov. = v, df = df, ...)
}

coeftest.survreg <- function(x, vcov. = NULL, df = Inf, ...)
{
  if(is.null(vcov.)) v <- vcov(x) else {
      if(is.function(vcov.)) v <- vcov.(x)
  	else v <- vcov.
  }
  if(length(x$coefficients) < NROW(x$var)) {
    x$coefficients <- c(x$coefficients, "Log(scale)" = log(x$scale))
  }
  coeftest.default(x, vcov. = v, df = df, ...)  
} 

coeftest.breakpointsfull <- function(x, vcov. = NULL, df = NULL, ..., save = FALSE)
{
  est <- coef(x, ...)
  if(is.null(df)) {
    df <- df.residual(x, ...)
    df <- as.vector(rep(df, rep(NCOL(est), length(df))))
  }  

  rnames <- as.vector(t(outer(rownames(est), colnames(est), paste)))
  est <- as.vector(t(est))
  
  se <- vcov(x, vcov. = vcov., ...)

  se <- as.vector(sapply(seq_along(se), function(x) sqrt(diag(se[[x]]))))
  tval <- est/se

  if(any(is.finite(df) && df > 0)) {
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
  attr(rval, "method") <- paste(mthd, "test of coefficients")
  ##  dQuote(class(x)[1]), "object", sQuote(deparse(substitute(x))))

  ## supplementary information for model summary
  attr(rval, "df") <- df
  attr(rval, "nobs") <- x$nobs
  attr(rval, "logLik") <- logLik(x, ...)
  if(save) attr(rval, "object") <- x

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

coef.coeftest <- function(object, ...) {
  object[, 1L, drop = TRUE]
}

df.residual.coeftest <- function(object, ...) {
  df <- attr(object, "df")
  if(df > 0) df else NULL
}

nobs.coeftest <- function(object, ...) {
  nobs <- attr(object, "nobs")
  if(nobs >= 0) nobs else NULL
}

logLik.coeftest <- function(object, ...) {
  attr(object, "logLik")
}

confint.coeftest <- function(object, parm = NULL, level = 0.95, ...)
{
  ## get estimates
  est <- object[, 1L]
  se <- object[, 2L]

  ## process level
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  
  ## get quantile from central limit theorem
  df <- attr(object, "df")
  if(is.null(df)) df <- 0
  fac <- if(is.finite(df) && df > 0) qt(a, df = df) else qnorm(a)

  ## set up confidence intervals
  ci <- cbind(est + fac[1] * se, est + fac[2] * se)
  colnames(ci) <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3L), "%")
  
  ## process parm
  if(is.null(parm)) parm <- seq_along(est)
  if(is.character(parm)) parm <- which(names(est) %in% parm)
  ci <- ci[parm, , drop = FALSE]
  return(ci)
} 
