reset <- function(formula,g, data=list())
{
	mf <- model.frame(formula, data=data);
	y <- model.response(mf);
	XX <- model.matrix(formula, data=data);
	n <- length(y);
	k <- length(XX[1,]);
	X <- XX[,1:g];
	Z <- XX[,(g+1):k];
	mod1 <- lm(y~X);
	mod2 <- lm(y~X+Z);
	beta1 <- coef(mod1);
	beta2 <- coef(mod2);
	res1 <- sum(resid(mod1)^2);
	res2 <- sum(resid(mod2)^2);
	df1 <- k-g;
	df2 <- n-k;
	pruefgroesse <- (res1 - res2) / res2;
	pruefgroesse <- (df2/df1)*pruefgroesse;
	names(pruefgroesse) <- "RESET";
	freiheit <- c(df1, df2);
	names(freiheit) <- c("df1","df2");
	RVAL <- list(statistic = pruefgroesse, 
			parameter = freiheit,
			method = "RESET-Test",
			p.value= 1-pf(pruefgroesse, df1, df2),
			data.name="formula");
	class(RVAL) <- "htest";
	return(RVAL);
};	


