hmctest <- function(formula, T, data=list())
{
	mf <- model.frame(formula, data=data);
	y <- model.response(mf);
	X <- model.matrix(formula, data=data);
	n <- length(y);
	k <- length(X[1,]);
	mod1 <- lm(y[1:T] ~ X[1:T,]);
	mod2 <- lm(y[(T+1):length(y)] ~ X[(T+1):length(y)]);
	resi1 <- resid(mod1);
	resi2 <- resid(mod2);
	sresi1 <- sum(resi1^2);
	sresi2 <- sum(resi2^2);
	hmc <- (sresi1/T)/(sresi2/n);
	names(hmc) <- "HMC";
	RVAL <- list(statistic = hmc, 
			method = "Harrison-McCabe-Test",
			p.value= NA,
			data.name="formula");
	class(RVAL) <- "htest";

	return(RVAL);
};	