gqtest <- function(formula, T, data=list())
{
	mf <- model.frame(formula, data=data);
	y <- model.response(mf);
	X <- model.matrix(formula, data=data);
	n <- length(y);
	k <- length(X[1,]);
	mod1 <- lm(y[1:T] ~ X[1:T,]);
	mod2 <- lm(y[(T+1):length(y)] ~ X[(T+1):length(y),]);
	resi1 <- resid(mod1);
	resi2 <- resid(mod2);
	sresi1 <- sum(resi1^2);
	sresi2 <- sum(resi2^2);
	gq <- (sresi2/(n-T-k))/(sresi1/(T-k));
	names(gq) <- "GQ";
	freiheit <- c(T-k, n-T-k);
	names(freiheit) <- c("df1", "df2");
	RVAL <- list(statistic = gq, 
			parameter = freiheit,
			method = "Goldfeld-Quandt-Test",
			p.value= 1-pf(gq, n-T-k, T-k),
			data.name="form1");
	class(RVAL) <- "htest";

	return(RVAL);
};	