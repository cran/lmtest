\name{reset}
\title{Reset-Test}
\usage{
rest(formula, g, data=list())
}
\alias{reset}
\arguments{
 \item{formula}{a symbolic describtion for the model to be tested}
 \item{g}{an integer indicating the number of regressors in model under H0}
 \item{data}{an optional data frame containing the variables in the model.
By default the variables are taken from the environment which 'reset' is
called from}
}
\description{
 \code{reset} performs the Ramsey's RESET-Test on any functional relation.
}
\details{
RESET-Test is popular means of diagnostic for correctness of functional
form. The basic assumption is that under the alternative, the model can be
written by the regression
\eqn{ y = X\beta + Z\gamma + u}{y=X * beta + Z * gamma}
\code{g} is the number of regressors in \code{X}. Ramsey's RESET is
equivalent  to the standard F-Test. 
The test statistic under \eqn{H_0} follows a F distribution with df1 and df2 degree
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
\author{Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de> and 
Frederik von Ameln <vonameln@amadeus.statistik.uni-dortmund.de>}

\seealso{\code{\link{lm}}}

\examples{
x <- c(1:30);
z <- x^2;
M <- cbind(x,z);
y <- M + rnorm(30,0,2);
formular <- y ~ M;
rain <- reset(formular,1);

\keyword{linear models}
