custest <- function(formula, data=list())
{
	mf <- model.frame(formula, data=data);
	y <- model.response(mf);
	X <- model.matrix(formula, data=data);
	k <- length(X[1,]);
	n <- length(y);
	mod <- lm(formula, data=data);
	resi <- resid(mod);
	sq <- sqrt(sum(resi^2)/(n - k));
	cu <- 0;
	cuplot <- c(0);
	for ( i in (k:(n-1)))
	{
		mod <- lm(y[1:i] ~ X[1:i,]-1);
		beta <- coef(mod);
		DUM <- solve(t(X[1:i,])%*%X[1:i,]);
		d1 <- y[i+1] - t(X[i+1,])%*%beta;
		d2 <- sqrt(1+t(X[i+1,])%*%DUM%*%X[i+1,]);
		u <- (1/sq)*(d1/d2);
		cu <- cu + u;
		cuplot <- c(cuplot, cu);
	}
	return(cuplot);
};