\name{raintest}
\title{Rainbow-Test}
\usage{
raintest(formula, data=list())
}
\alias{raintest}
\arguments{
 \item{formula}{a symbolic describtion for the model to be tested}
 \item{data}{an optional data frame containing the variables in the model.
By default the variables are taken from the environment which 'raintest' is
called from}
}
\description{
 \code{raintest} performs the Rainbow-Tests on linearity. 
}
\details{
The basic idea of the Rainbow-Test is that even if the true relationship is
nonlinear, over a subsample of data given a good linear fit can be achieved.
The null hypothesis is rejected whenever the overall fit is significantly
inferious to the fit of the subsample. The test statistic under \eqn{H_0}
follows a F distribution with df1 and df2 degree of freedom. This
particular procedure compares a subsample consisting of all data points
without the upper and lower quartile (of time index), thus 50% of all data points, with all
data points. If the true relationship is concave or convex, the null hypothesis
should be rejected. }
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
x <- c(1:30)
y <- x^2 + rnorm(30,0,2)
formular <- y ~ x
rain <- raintest(formular)
qf(0.95, rain$parameter[1], rain$parameter[2])
}
\keyword{htest}
