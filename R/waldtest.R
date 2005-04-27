waldtest <- function(object, ...) {
  UseMethod("waldtest")
}

waldtest.formula <- function(object, ...) {
  waldtest.lm(object, ...)
}

waldtest.lm <- function(object, ..., vcov = NULL, test = c("F", "Chisq"), data = list())
{
  ## convenience functions:
  ## 1. extracts term labels
  tlab <- function(x) attr(terms(x), "term.labels")
  ## 2. compute an updated "lm" object
  lmUpdate <- function(fm, update) {
    ## if `update' is numeric or character, then assume that the 
    ## corresponding variables (i.e., terms) are redundant (i.e., should be omitted)
    if(is.numeric(update)) {
      ## sanity checking of numeric update specification
      if(any(update < 1)) {
        warning("for numeric model specifications all values have to be >=1")
	update <- abs(update)[abs(update) > 0]
      }
      if(any(update > length(tlab(fm)))) {
        warning(paste("more terms specified than existent in the model:",
	        paste(as.character(update[update > length(tlab(fm))]), collapse = ", ")))
	update <- update[update <= length(tlab(fm))]
      }
      ## finally turn numeric into character update specification
      update <- tlab(fm)[update]
    }
    if(is.character(update)) {
      ## sanity checking of character update specification
      if(!all(update %in% tlab(fm))) {
        warning(paste("terms specified that are not in the model:",
	        paste(dQuote(update[!(update %in% tlab(fm))]), collapse = ", ")))
        update <- update[update %in% tlab(fm)]
      }
      if(length(update) < 1) stop("empty model specification")  
      ## finally turn character into formula update specification       
      update <- as.formula(paste(". ~ . -", paste(update, collapse = " - ")))
    }
    if(inherits(update, "formula")) update <- update(fm, update)
    if(!inherits(update, "lm")) stop("invalid update model specification")
    return(update)
  }
  ## 3. compare two "lm" objects
  lmCompare <- function(fm, fm.up, vfun = NULL) {
    q <- length(coef(fm)) - length(coef(fm.up))

    if(q > 0) {
      fm0 <- fm.up
      fm1 <- fm
    } else {
      fm0 <- fm
      fm1 <- fm.up
    }
    k <- length(coef(fm1))
    n <- length(fm1$residuals)

    ## determine omitted variables
    if(!all(tlab(fm0) %in% tlab(fm1))) stop("models are not nested")
    ovar <- which(!(names(coef(fm1)) %in% names(coef(fm0))))

    ## get covariance matrix estimate
    vc <- if(is.null(vfun)) vcov(fm1)
          else if(is.function(vfun)) vfun(fm1)
	  else vfun
  
    ## compute Chisq statistic
    stat <- t(coef(fm1)[ovar]) %*% solve(vc[ovar,ovar]) %*% coef(fm1)[ovar]
    return(c(q, stat))
  }

  ## check whether `object' is already a fitted model
  if(!inherits(object, "lm")) object <- lm(object, data = data)

  ## recursively fit all objects (if necessary)
  objects <- list(object, ...)
  nmodels <- length(objects)
  if(nmodels < 2) {
    objects <- c(objects, . ~ 1)
    nmodels <- 2
  }
  for(i in 2:nmodels) objects[[i]] <- lmUpdate(objects[[i-1]], objects[[i]])

  ## check responses
  responses <- as.character(lapply(objects, function(x) deparse(terms(x)[[2]])))
  sameresp <- responses == responses[1]
  if(!all(sameresp)) {
    objects <- objects[sameresp]
    warning("models with response ", deparse(responses[!sameresp]),
	    " removed because response differs from ", "model 1")
  }
  ## check sample sizes
  ns <- sapply(objects, function(x) length(x$residuals))
  if(any(ns != ns[1])) stop("models were not all fitted to the same size of dataset")
  ## check vcov
  if(nmodels > 2 && !is.null(vcov) && !is.function(vcov))
    stop("to compare more than 2 models `vcov' needs to be a function")
  
  ## setup ANOVA matrix
  test <- match.arg(test)
  rval <- matrix(rep(NA, 4 * nmodels), ncol = 4)
  colnames(rval) <- c("Res.Df", "Df", test, paste("Pr(>", test, ")", sep = ""))
  rownames(rval) <- 1:nmodels
  rval[,1] <- as.numeric(sapply(objects, df.residual))
  for(i in 2:nmodels) rval[i, 2:3] <- lmCompare(objects[[i-1]], objects[[i]], vfun = vcov)
  if(test == "Chisq") {
    rval[,4] <- pchisq(rval[,3], abs(rval[,2]), lower.tail = FALSE)
  } else {
    df <- rval[,1]
    for(i in 2:nmodels) if(rval[i,2] > 0) df[i] <- df[i-1]
    rval[,3] <- rval[,3]/abs(rval[,2])
    rval[,4] <- pf(rval[,3], abs(rval[,2]), df, lower.tail = FALSE)
  }

  variables <- lapply(objects, function(x) paste(deparse(formula(x)), collapse="\n") )
  title <- "Wald test\n"
  topnote <- paste("Model ", format(1:nmodels),": ", variables, sep="", collapse="\n")

  structure(as.data.frame(rval), heading = c(title, topnote),
	    class = c("anova", "data.frame"))
}
