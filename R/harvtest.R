harvtest <- function(formula, data=list())
{
	mf <- model.frame(formula, data=data);
	y <- model.response(mf);
	X <- model.matrix(formula, data=data);
	k <- length(X[1,]);
	n <- length(y);

	mod <- lm(formula, data=data);
	resi <- resid(mod);
	cu <- 0;
	cuplot <- c(0);
	rekres <- c(0);
	sq <- sqrt(sum(resi^2)/(n - k));
	
	for ( i in ((k+1):(n-1)))
	{
		mod <- lm(y[1:i] ~ X[1:i,]-1);
		beta <- coef(mod);
		DUM <- solve(t(X[1:i,])%*%X[1:i,]);
		d1 <- y[i+1] - t(X[i+1,])%*%beta;
		d2 <- sqrt(1+t(X[i+1,])%*%DUM%*%X[i+1,]);
		u <- (1/sq)*(d1/d2);
		cu <- cu + u;
		rekres <- c(rekres,u);
		cuplot <- c(cuplot, cu);
	}
	
	reksum <- sum((rekres - mean(rekres))^2);
	harv <- abs((cu/sqrt(n -k))/sqrt(reksum/(n-k-1)));
	names(harv) <- "HC";
	freiheit <- n-k-1;
	names(freiheit) <- "df";
	RVAL <- list(statistic = harv, 
			parameter = freiheit,
			method = "Harvey-Collier-Test",
			p.value= 1-pt(harv, n-k-1),
			data.name="formula");
	class(RVAL) <- "htest";

	return(RVAL);
};