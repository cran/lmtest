\name{harvtest}
\title{Harvey-Collier-Test}
\usage{
harvtest(formula, data=list())
}
\alias{harvtest}
\arguments{
 \item{formula}{a symbolic describtion for the model to be tested}
 \item{data}{an optional data frame containing the variables in the model.
By default the variables are taken from the environment which 'harvtest' is
called from}
}
\description{
 \code{harvtest} performs the Harvey-Collier-Test for linearity.
}
\details{
The Harvey-Collier-Tests examines if the true relationship is linear versus a
convex or concave relationship under the alternative. The test statistic under
null hypothesis follows a t distribution with df degree of freedom.
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
x <- c(1:30);
y <- x^2 + rnorm(30,0,2);
formular <- y ~ x;
harv <- harvtest(formular);
qt(0.95, harv$parameter[1]);
}
\keyword{linear models}ÿ
