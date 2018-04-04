\name{moneydemand}
\alias{moneydemand}
\docType{data}
\encoding{latin1}
\title{Money Demand}
\usage{data(moneydemand)}
\description{
Money Demand Data.
}
\format{
 A multivariate yearly time series from 1879 to 1974 with variables
 \describe{
    \item{logM}{logarithm of quantity of money,}
    \item{logYp}{logarithm of real permanent income,}
    \item{Rs}{short term interest rate,}
    \item{Rm}{rate of return on money,}
    \item{Rl}{not documented in the sources,}
    \item{logSpp}{logarithm of an operational measure of the variability of
the rate of price changes.}
 }
}
\source{
The data was originally studied by Allen (1982), the data set is given in
Krämer and Sonnberger (1986). Below we replicate a few examples 
from the book.
Some of these results differ more or less seriously and are sometimes
parameterized differently.

}
\references{

S.D. Allen (1982),
Klein's Price Variability Terms in the U.S. Demand for Money.
\emph{Journal of Money, Credit and Banking} \bold{14}, 525--530

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica


}
\examples{
data(moneydemand)
moneydemand <- window(moneydemand, start=1880, end=1972)

## page 125, fit Allen OLS model (and Durbin-Watson test),
## last line in Table 6.1

modelAllen <- logM ~ logYp + Rs + Rl + logSpp
lm(modelAllen, data = moneydemand)
dwtest(modelAllen, data = moneydemand)

## page 127, fit test statistics in Table 6.1 c)
################################################

## Breusch-Pagan
bptest(modelAllen, studentize = FALSE, data = moneydemand)
bptest(modelAllen, studentize = TRUE, data = moneydemand)

## RESET
reset(modelAllen, data = moneydemand)
reset(modelAllen, power = 2, type = "regressor", data = moneydemand)
reset(modelAllen, type = "princomp", data = moneydemand)

## Harvey-Collier tests (up to sign of the test statistic)
harvtest(modelAllen, order.by = ~logYp, data = moneydemand)
harvtest(modelAllen, order.by = ~Rs, data = moneydemand)
harvtest(modelAllen, order.by = ~Rl, data = moneydemand)
harvtest(modelAllen, order.by = ~logSpp, data = moneydemand)

## Rainbow test
raintest(modelAllen, order.by = "mahalanobis", data = moneydemand)


if(require(strucchange, quietly = TRUE)) {
## Chow (1913)
sctest(modelAllen, point=c(1913,1), data = moneydemand, type = "Chow") }

if(require(strucchange, quietly = TRUE)) {
## Fluctuation
sctest(modelAllen, type = "fluctuation", rescale = FALSE, data = moneydemand)}
}
\keyword{datasets}
