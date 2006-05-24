lrtest <- function(object, ...) {
  UseMethod("lrtest")
}

lrtest.formula <- function(object, ..., data = list()) {
  object <- if(length(data) < 1) eval(call("lm", formula = as.formula(deparse(substitute(object))),
    environment(object)))
  else eval(call("lm", formula = as.formula(deparse(substitute(object))),
    data = as.name(deparse(substitute(data))), environment(data)))
  lrtest.default(object, ...)
}

lrtest.default <- function(object, ..., name = NULL)
{
  ## methods needed:
  ## - logLik()

  ## - terms()  -> for updating only
  ## - update()

  ## - formula() -> for determining `name' (can be user-supplied)
  ## - residuals() -> only for determining number of observations

  ## model class
  cls <- class(object)[1]

  ## convenience functions:
  ## 1. extract number of observations
  nobs <- function(x) NROW(residuals(x))
  ## 2. extracts term labels
  tlab <- function(x) attr(terms(x), "term.labels")
  ## 3. extracts model name
  if(is.null(name)) name <- function(x) paste(deparse(formula(x)), collapse="\n")
  ## 3. compute an updated model object
  modelUpdate <- function(fm, update) {
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
    if(!inherits(update, cls)) stop(paste("original model was of class \"", cls,
      "\", updated model is of class \"", class(update)[1], "\"", sep = ""))
    return(update)
  }

  ## recursively fit all objects (if necessary)
  objects <- list(object, ...)
  nmodels <- length(objects)
  if(nmodels < 2) {
    objects <- c(objects, . ~ 1)
    nmodels <- 2
  }
  
  # remember which models are already fitted and which are described
  # by an update mechanism
  no.update <- sapply(objects, function(obj) inherits(obj, cls))
  
  ## updating
  for(i in 2:nmodels) objects[[i]] <- modelUpdate(objects[[i-1]], objects[[i]])

  ## check responses
  responses <- as.character(lapply(objects, function(x) deparse(terms(x)[[2]])))
  sameresp <- responses == responses[1]
  if(!all(sameresp)) {
    objects <- objects[sameresp]
    warning("models with response ", deparse(responses[!sameresp]),
	    " removed because response differs from ", "model 1")
  }

  ## check sample sizes
  ns <- sapply(objects, nobs)
  if(any(ns != ns[1])) {
    for(i in 2:nmodels) {
      if(ns[1] != ns[i]) {
        if(no.update[i]) stop("models were not all fitted to the same size of dataset")
	  else {
	    commonobs <- row.names(model.frame(objects[[i]])) %in% row.names(model.frame(objects[[i-1]]))
	    objects[[i]] <- eval(substitute(update(objects[[i]], subset = commonobs),
	      list(commonobs = commonobs)))
	    if(nobs(objects[[i]]) != ns[1]) stop("models could not be fitted to the same size of dataset")
	  }
      }
    }
  }

  ## setup ANOVA matrix
  rval <- matrix(rep(NA, 5 * nmodels), ncol = 5)
  colnames(rval) <- c("#Df", "LogLik", "Df", "Chisq", "Pr(>Chisq)")
  rownames(rval) <- 1:nmodels
  
  logL <- lapply(objects, logLik)
  rval[,1] <- as.numeric(sapply(logL, function(x) attr(x, "df")))  
  rval[,2] <- sapply(logL, as.numeric)
  rval[2:nmodels, 3] <- rval[2:nmodels, 1] - rval[1:(nmodels-1), 1]
  rval[2:nmodels, 4] <- 2 * abs(rval[2:nmodels, 2] - rval[1:(nmodels-1), 2])
  rval[,5] <- pchisq(rval[,4], round(abs(rval[,3])), lower.tail = FALSE)

  variables <- lapply(objects, name)
  title <- "Likelihood ratio test\n"
  topnote <- paste("Model ", format(1:nmodels),": ", variables, sep="", collapse="\n")

  structure(as.data.frame(rval), heading = c(title, topnote),
	    class = c("anova", "data.frame"))
}
