\name{harvtest}
\alias{harvtest}
\encoding{latin1}
\title{Harvey-Collier Test}
\description{Harvey-Collier test for linearity.}

\usage{
harvtest(formula, order.by = NULL, data = list())
}

\arguments{
  \item{formula}{a symbolic description for the model to be tested
    (or a fitted \code{"lm"} object).}
 \item{order.by}{Either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series).}
 \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which
    \code{harvtest} is called from.}
}

\details{
The Harvey-Collier test performs a t-test (with \code{parameter} degrees of
freedom) on the recursive residuals. If the true relationship is not linear but
convex or concave the mean of the recursive residuals should differ
from 0 significantly.

Examples can not only be found on this page, but also on the help pages of the
data sets \code{\link{bondyield}}, \code{\link{currencysubstitution}},
\code{\link{growthofmoney}}, \code{\link{moneydemand}},
 \code{\link{unemployment}},
\code{\link{wages}}.
}

\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{parameter}{degrees of freedom.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}

\references{

A. Harvey & P. Collier (1977),
Testing for Functional Misspecification in Regression Analysis.
\emph{Journal of Econometrics} \bold{6}, 103--119

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica


}

\seealso{\code{\link{lm}}}

\examples{
# generate a regressor and dependent variable
x <- 1:50
y1 <- 1 + x + rnorm(50)
y2 <- y1 + 0.3*x^2

## perform Harvey-Collier test
harv <- harvtest(y1 ~ x)
harv
## calculate critical value vor 0.05 level
qt(0.95, harv$parameter)
harvtest(y2 ~ x)
}
\keyword{htest}
