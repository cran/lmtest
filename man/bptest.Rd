\name{bptest}
\title{Breusch-Pagan-Test}
\usage{
bptest(formula, data=list())
}
\alias{bptest}
\arguments{
 \item{formula}{a symbolic describtion for the model to be tested}
 \item{data}{an optional data frame containing the variables in the model.
By default the variables are taken from the environment which 'bptest' is
called from}
}
\description{
 \code{bptest} performs the Breusch-Pagan-Test against heteroskedasticity. 
}
\details{
The test statistic under H0 follows a chisq distribution with df degree
of freedom.
}
\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Kraemer, W., Sonnberger, H. (1986): The linear regression model
under test}
\author{Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de>}

\seealso{\code{\link{lm}}}

\examples{
# generate a regressor;
x <- sort(runif(50,-10,10));
# generate heteroskedastic disturbances
err <- c(rnorm(25,0,1),rnorm(25,0,10));
# generate a linear relationship
y <- 2*x + err;
# perform Breusch-Pagan
bp <- bptest(y~x-1);
# calculate quantil
qchisq(0.95, bp$parameter[1]);
}
\keyword{htest}