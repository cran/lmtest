\name{bgtest}
\alias{bgtest}
\title{Breusch-Godfrey Test}
\description{
 \code{bgtest} performs the Breusch-Godfrey test for higher-order
 serial correlation.
}
\usage{
bgtest(formula, order = 1, order.by = NULL, type = c("Chisq", "F"), data = list())
}
\arguments{
  \item{formula}{a symbolic description for the model to be tested
    (or a fitted \code{"lm"} object).}
  \item{order}{integer. maximal order of serial correlation to be tested.}
 \item{order.by}{Either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series).}
  \item{type}{the type of test statistic to be returned. Either
    \code{"Chisq"} for the Chi-squared test statistic or
    \code{"F"} for the F test statistic.}
  \item{data}{an optional data frame containing the variables in the
    model. By default the variables are taken from the environment
    which \code{bgtest} is called from.}
  }
\details{
  Under \eqn{H_0} the test statistic is asymptotically Chi-squared with
  degrees of freedom as given in \code{parameter}.
  If \code{type} is set to \code{"F"} the function returns
  the exact F statistic which, under \eqn{H_0}, follows an \eqn{F}
  distribution with degrees of freedom as given in \code{parameter}.

  The starting values for the lagged residuals in the supplementary
  regression are chosen to be 0.
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
  Johnston, J. (1984): \emph{Econometric Methods}, Third Edition, McGraw Hill
  Inc.

  Godfrey, L.G. (1978): `Testing Against General Autoregressive and
  Moving Average Error Models when the Regressors Include Lagged
  Dependent Variables', \emph{Econometrica}, 46, 1293-1302.

  Breusch, T.S. (1979): `Testing for Autocorrelation in Dynamic Linear
  Models', \emph{Australian Economic Papers}, 17, 334-355.
}
\author{David Mitchell <david.mitchell@dotars.gov.au>, Achim Zeileis}
\seealso{\code{\link{dwtest}}}
\examples{

     ## Generate a stationary and an AR(1) series
     x <- rep(c(1, -1), 50)

     y1 <- 1 + x + rnorm(100)

     ## Perform Breusch-Godfrey test for first-order serial correlation:
     bgtest(y1 ~ x)
     ## or for fourth-order serial correlation
     bgtest(y1 ~ x, order = 4)
     ## Compare with Durbin-Watson test results:
     dwtest(y1 ~ x)

     y2 <- filter(y1, 0.5, method = "recursive")
     bgtest(y2 ~ x)

}
\keyword{htest}
