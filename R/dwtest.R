dwtest <- function(formula, data=list())
{
	mod <- lm(formula, data=data);
	resi <- resid(mod);
	dw <- sum((resi[2:length(resi)] - resi[1:length(resi)-1])^2)/sum(resi^2);
	names(dw) <- "DW";
	RVAL <- list(statistic = dw, 
			method = "Durbin-Watson-Test",
			p.value= 0,
			data.name=" ");
	class(RVAL) <- "htest";

	return(RVAL);
};
