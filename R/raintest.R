raintest <- function(formula, data=list())
{
	mf <- model.frame(formula, data=data);
	y <- model.response(mf);
	X <- model.matrix(formula, data=data);
	n <- length(y);
	k <- length(X[1,]);
	g1 <- ceiling(quantile(c(1:n), probs=0.25));
	g2 <- floor(quantile(c(1:n), probs=0.75));
	n1 <- g2 - g1;
	y1 <- y[g1:g2];
	X1 <- X[g1:g2,];
	mod <- lm(formula, data=data);
	mod1 <- lm(y1 ~X1);
	sresi <- sum(resid(mod)^2);
	sresi1 <- sum(resid(mod1)^2);
	rain <- ((sresi - sresi1)/(n-n1))/(sresi1/(n1-k));	
	names(rain) <- "Rain";
	freiheit <- c((n-n1),(n1-k));
	names(freiheit) <- c("df1","df2");
	RVAL <- list(statistic = rain, 
			parameter = freiheit,
			method = "Rainbow-Test",
			p.value= 1-pf(rain, n-n1, n1-k),
			data.name="formula");
	class(RVAL) <- "htest";

	return(RVAL);
};	

