\name{currencysubstitution}
\alias{currencysubstitution}
\docType{data}
\encoding{latin1}
\title{ Currency Substitution }
\usage{data(currencysubstitution)}
\description{
Currency Substitution Data.
}
\format{
 A multivariate quarterly time series from 1960(4) to 1975(4) with variables
 \describe{
    \item{logCUS}{logarithm of the ratio of Canadian holdings of Canadian
dollar balances and Canadian holdings of U.S. dollar balances,}
    \item{Iu}{yield on U.S. Treasury bills,}
    \item{Ic}{yield on Canadian Treasury bills,}
    \item{logY}{logarithm of Canadian real gross national product.}
  }
}
\source{
The data was originally studied by Miles (1978), the data set is given in
Krämer and Sonnberger (1986). Below we replicate a few examples 
from their book.
Some of these results differ more or less seriously and are sometimes
parameterized differently.
}
\references{

M. Miles (1978),
Currency Substitution, Flexible Exchange Rates, and Monetary Independence.
\emph{American Economic Review} \bold{68}, 428--436

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica


}
\examples{
data(currencysubstitution)

## page 130, fit Miles OLS model and Bordo-Choudri OLS model
## third and last line in Table 6.3

modelMiles <- logCUS ~ log((1+Iu)/(1+Ic))
lm(modelMiles, data=currencysubstitution)
dwtest(modelMiles, data=currencysubstitution)

modelBordoChoudri <- logCUS ~ I(Iu-Ic) + Ic + logY
lm(modelBordoChoudri, data=currencysubstitution)
dwtest(modelBordoChoudri, data=currencysubstitution)


## page 131, fit test statistics in Table 6.4 b)
################################################

if(require(strucchange, quietly = TRUE)) {
## Fluctuation test
sctest(modelMiles, type="fluctuation", data=currencysubstitution,
 rescale=FALSE) }

## RESET
reset(modelMiles, data=currencysubstitution)
reset(modelMiles, power=2, type="regressor", data=currencysubstitution)
reset(modelMiles, type="princomp", data=currencysubstitution)

## Harvey-Collier
harvtest(modelMiles, order.by = ~log((1+Iu)/(1+Ic)), data=currencysubstitution)

## Rainbow
raintest(modelMiles, order.by = "mahalanobis", data=currencysubstitution)


## page 132, fit test statistics in Table 6.4 d)
################################################

if(require(strucchange, quietly = TRUE)) {
## Chow 1970(2)
sctest(modelBordoChoudri, point=c(1970,2), data=currencysubstitution,
 type="Chow") }

## Breusch-Pagan
bptest(modelBordoChoudri, data=currencysubstitution, studentize=FALSE)
bptest(modelBordoChoudri, data=currencysubstitution)

## RESET
reset(modelBordoChoudri, data=currencysubstitution)
reset(modelBordoChoudri, power=2, type="regressor", data=currencysubstitution)
reset(modelBordoChoudri, type="princomp", data=currencysubstitution)

## Harvey-Collier
harvtest(modelBordoChoudri, order.by = ~ I(Iu-Ic), data=currencysubstitution)
harvtest(modelBordoChoudri, order.by = ~ Ic, data=currencysubstitution)
harvtest(modelBordoChoudri, order.by = ~ logY, data=currencysubstitution)

## Rainbow
raintest(modelBordoChoudri, order.by = "mahalanobis", data=currencysubstitution)
}
\keyword{datasets}
