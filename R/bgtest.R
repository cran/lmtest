bgtest <- function(formula, order = 1, order.by = NULL, type = c("Chisq", "F"), data = list())
{
  dname <- paste(deparse(substitute(formula)))

  if(!inherits(formula, "formula")) {
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
