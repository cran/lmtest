\name{hmctest}
\alias{hmctest}
\encoding{latin1}
\title{Harrison-McCabe test}
\description{Harrison-McCabe test for heteroskedasticity.}
\usage{
hmctest(formula, point = 0.5, order.by = NULL, simulate.p = TRUE, nsim = 1000,
  plot = FALSE, data = list())
}

\arguments{
  \item{formula}{a symbolic description for the model to be tested
    (or a fitted \code{"lm"} object).}
 \item{point}{numeric. If \code{point} is smaller than 1 it is
     interpreted as percentages of data, i.e. \code{n*point} is
     taken to be the (potential) breakpoint in the variances, if
     \code{n} is the number of observations in the model. If \code{point}
     is greater than 1 it is interpreted to be the index of the breakpoint.}
 \item{order.by}{Either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series).}
  \item{simulate.p}{logical. If \code{TRUE} a p value will be assessed
    by simulation, otherwise the p value is \code{NA}.}
  \item{nsim}{integer. Determines how many runs are used to simulate the
    p value.}
  \item{plot}{logical. If \code{TRUE} the test statistic for all possible
    breakpoints is plotted.}
  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which \code{hmctest}
    is called from.}
}

\details{The Harrison-McCabe test statistic is the fraction of the residual sum
 of
squares that relates to the fraction of the data before the breakpoint. Under
 \eqn{H_0}
the test statistic should be close to the size of this fraction, e.g. in the
 default
case close to 0.5. The null hypothesis is reject if the statistic is too small.

Examples can not only be found on this page, but also on the help pages of the
data sets \code{\link{bondyield}}, \code{\link{currencysubstitution}},
\code{\link{growthofmoney}}, \code{\link{moneydemand}},
 \code{\link{unemployment}},
\code{\link{wages}}.}

\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the simulated p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}

\references{

M.J. Harrison & B.P.M McCabe (1979),
A Test for Heteroscedasticity based on Ordinary Least Squares Residuals.
\emph{Journal of the American Statistical Association} \bold{74}, 494--499

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica


}

\seealso{\code{\link{lm}}}

\examples{
## generate a regressor
x <- rep(c(-1,1), 50)
## generate heteroskedastic and homoskedastic disturbances
err1 <- c(rnorm(50, sd=1), rnorm(50, sd=2))
err2 <- rnorm(100)
## generate a linear relationship
y1 <- 1 + x + err1
y2 <- 1 + x + err2
## perform Harrison-McCabe test
hmctest(y1 ~ x)
hmctest(y2 ~ x)
}
\keyword{htest}
