\name{wages}
\alias{wages}
\docType{data}
\encoding{latin1}
\title{Wages}
\usage{data(wages)}
\description{
Wages Data.
}
\format{
 A multivariate yearly time series from 1960 to 1979 with variables
 \describe{
    \item{w}{wages,}
    \item{CPI}{consumer price index,}
    \item{u}{unemployment,}
    \item{mw}{minimum wage.}
  }

}
\source{The data was originally studied by Nicols (1983), the data set is given
in Krämer and Sonnberger (1986). Below we replicate a few examples from their
 book.
Some of these results differ more or less seriously and are sometimes
parameterized differently.
}
\references{

D.A. Nicols (1983),
Macroeconomic Determinants of Wage Adjustments in White Collar Occupations.
\emph{Review of Economics and Statistics} \bold{65}, 203--213

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica


}
\examples{
data(wages)

## data transformation to include lagged series
mywages <- cbind(wages, lag(wages[,2], k = -1), lag(wages[,2], k = -2))
colnames(mywages) <- c(colnames(wages), "CPI2", "CPI3")
mywages <- window(mywages, start=1962, end=1979)

## page 142, fit Nichols OLS model
## equation (6.10)

modelNichols <- w ~ CPI + CPI2 + CPI3 + u + mw
lm(modelNichols, data = mywages)

## page 143, fit test statistics in table 6.11
##############################################

if(require(strucchange, quietly = TRUE)) {
## Chow 1972
sctest(modelNichols, point=c(1971,1), data=mywages, type="Chow") }

## Breusch-Pagan
bptest(modelNichols, data=mywages, studentize=FALSE)
bptest(modelNichols, data=mywages)

## RESET (a)-(b)
reset(modelNichols, data=mywages)
reset(modelNichols, power=2, type="regressor", data=mywages)

## Harvey-Collier
harvtest(modelNichols, order.by = ~ CPI, data=mywages)
harvtest(modelNichols, order.by = ~ CPI2, data=mywages)
harvtest(modelNichols, order.by = ~ CPI3, data=mywages)
harvtest(modelNichols, order.by = ~ u, data=mywages)

## Rainbow
raintest(modelNichols, order.by = "mahalanobis", data=mywages)
}
\keyword{datasets}
