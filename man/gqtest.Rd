\name{gqtest}
\alias{gqtest}
\encoding{latin1}
\title{Goldfeld-Quandt Test}
\description{Goldfeld-Quandt test against heteroskedasticity.}
\usage{
gqtest(formula, point = 0.5, fraction = 0,
  alternative = c("greater", "two.sided", "less"),
  order.by = NULL, data = list())
}
\arguments{
  \item{formula}{a symbolic description for the model to be tested
    (or a fitted \code{"lm"} object).}
 \item{point}{numerical. If \code{point} is smaller than 1 it is
     interpreted as percentages of data, i.e. \code{n*point} is
     taken to be the (potential) breakpoint in the variances, if
     \code{n} is the number of observations in the model. If \code{point}
     is greater than 1 it is interpreted to be the index of the breakpoint.}
 \item{fraction}{numerical. The number of central observations to be omitted.
    If \code{fraction} is smaller than 1, it is chosen to be \code{fraction*n}
    if \code{n} is the number of observations in the model.}
 \item{alternative}{a character string specifying the alternative hypothesis.
    The default is to test for increasing variances.}
 \item{order.by}{Either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series).}
 \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which \code{gqtest}
    is called from.}
}
\details{The Goldfeld-Quandt test compares the variances of two submodels
divided by a specified breakpoint and rejects if the variances differ.

Under \eqn{H_0} the test statistic of the Goldfeld-Quandt test follows an F
distribution with the degrees of freedom as given in \code{parameter}.

Examples can not only be found on this page, but also on the help pages of the
data sets \code{\link{bondyield}}, \code{\link{currencysubstitution}},
\code{\link{growthofmoney}}, \code{\link{moneydemand}},
 \code{\link{unemployment}}, \code{\link{wages}}.
}
\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{parameter}{degrees of freedom.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{p.value}{the p-value of the test.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{

S.M. Goldfeld & R.E. Quandt (1965),
Some Tests for Homoskedasticity.
\emph{Journal of the American Statistical Association} \bold{60}, 539--547

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
## perform Goldfeld-Quandt test
gqtest(y1 ~ x)
gqtest(y2 ~ x)
}
\keyword{htest}
