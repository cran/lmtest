bptest <- function(formula, data=list())
{
	mod <- lm(formula, data=data);
	X <- model.matrix(formula, data=data);
	Z <- X;
	n <- length(X[,1]);
	k <- length(X[1,]);
	resi <- resid(mod);
	s2 <- (t(resi)%*%resi)/(n-k); 
	T <- c(1:n)
	f <- c(resi[T]^2/s2 -1);
	bp <- 1/2*t(f)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%f;
	names(bp) <- "BP";
	freiheit <- c(k);
	names(freiheit) <- "df";
	RVAL <- list(statistic = bp, 
			parameter = freiheit,
			method = "Breusch-Pagan-Test",
			p.value= 1-pchisq(bp,k),
			data.name="form1");
	class(RVAL) <- "htest";

	return(RVAL);
};

