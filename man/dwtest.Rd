\name{dwtest}
\title{Durbin-Watson-Test}
\usage{
raintest(formula, data=list())
}
\alias{dwtest}
\arguments{
 \item{formula}{a symbolic describtion for the model to be tested}
 \item{data}{an optional data frame containing the variables in the model.
By default the variables are taken from the environment which 'dwtest' is
called from}
}
\description{
 \code{dwtest} performs the Durbin-Watson-Test on autocorrelation of
disturbances.
}
\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Kraemer, W., Sonnberger, H. (1986): The linear regression model
under test}
\author{Torsten Hothorn <hothorn@amadeus.statistik.uni-dortmund.de>}

\seealso{\code{\link{lm}}}

\examples{
x <- c(1:30);
ut <- c(1:30);
# creating a ARMA(1) process
for(i in (2:30))
{
	ut[i] <- 1 * ut[i-1] + rnorm(1,0,2);
}
y <- x + ut;
dw <- dwtest(y ~ x);
}
\keyword{linear models}