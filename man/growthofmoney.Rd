\name{growthofmoney}
\alias{growthofmoney}
\docType{data}
\encoding{latin1}
\title{Growth of Money Supply}
\usage{data(growthofmoney)}
\description{
Growth of Money Supply Data.
}
\format{
 A multivariate quarterly time series from 1970(2) to 1974(4) with variables
 \describe{
    \item{TG1.TG0}{difference of current and preceding target for the growth
 rate of the
money supply,}
    \item{AG0.TG0}{difference of actual growth rate and target growth
rate for the preceding period.}
 }
}
\source{
The data was originally studied by Hetzel (1981), the data set is given in
Krämer and Sonnberger (1986). Below we replicate a few examples from their
 book.
Some of these results differ more or less seriously and are sometimes
parameterized differently.
}
\references{

R.L. Hetzel (1981),
The Federal Reserve System and Control of the Money Supply in the 1970's.
\emph{Journal of Money, Credit and Banking} \bold{13}, 31--43

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica


}
\examples{
data(growthofmoney)

## page 137, fit Hetzel OLS model
## first/second line in Table 6.7

modelHetzel <- TG1.TG0 ~ AG0.TG0
lm(modelHetzel, data=growthofmoney)
dwtest(modelHetzel, data=growthofmoney)


## page 135, fit test statistics in Table 6.8
#############################################

if(require(strucchange, quietly = TRUE)) {
## Chow 1974(1)
sctest(modelHetzel, point=c(1973,4), data=growthofmoney, type="Chow") }

## RESET
reset(modelHetzel, data=growthofmoney)
reset(modelHetzel, power=2, type="regressor", data=growthofmoney)
reset(modelHetzel, type="princomp", data=growthofmoney)

## Harvey-Collier
harvtest(modelHetzel, order.by= ~ AG0.TG0, data=growthofmoney)

## Rainbow
raintest(modelHetzel, order.by = "mahalanobis", data=growthofmoney)


## Identification of outliers
#############################

## Figure 6.1
plot(modelHetzel, data=growthofmoney)
abline(v=0)
abline(h=0)
abline(coef(lm(modelHetzel, data=growthofmoney)), col=2)

## Table 6.7, last line
growthofmoney2 <- as.data.frame(growthofmoney[-c(5:6),])
lm(modelHetzel, data=growthofmoney2)
dwtest(modelHetzel, data=growthofmoney2)
}
\keyword{datasets}
