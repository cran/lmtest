\name{unemployment}
\alias{unemployment}
\docType{data}
\encoding{latin1}
\title{Unemployment Data}
\usage{data(unemployment)}
\description{
Unemployment Data.
}
\format{
 A multivariate yearly time series from 1890 to 1979 with variables
 \describe{
    \item{UN}{unemployment rate,}
    \item{m}{broad money supply,}
    \item{p}{implicit deflator of Gross National Product,}
    \item{G}{real purchases of goods and services,}
    \item{x}{real exports.}
  }

}
\source{
The data was originally studied by Rea (1983), the data set is given in
Krämer and Sonnberger (1986). Below we replicate a few examples 
from their book. Some of these results differ more or less seriously and are
 sometimes
parameterized differently.
}
\references{

J.D. Rea (1983),
The Explanatory Power of Alternative Theories of Inflation and Unemployment,
1895-1979.
\emph{Review of Economics and Statistics} \bold{65}, 183--195

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica


}
\examples{
data(unemployment)

## data transformation
myunemployment <- window(unemployment, start=1895, end=1956)
time <- 6:67

## page 144, fit Rea OLS model
## last line in Table 6.12

modelRea <- UN ~ log(m/p) + log(G) + log(x) + time
lm(modelRea, data = myunemployment)
## coefficients of logged variables differ by factor 100

## page 143, fit test statistics in table 6.11
##############################################

if(require(strucchange, quietly = TRUE)) {
## Chow 1941
sctest(modelRea, point=c(1940,1), data=myunemployment, type="Chow") }

## Breusch-Pagan
bptest(modelRea, data=myunemployment, studentize=FALSE)
bptest(modelRea, data=myunemployment)

## RESET (a)-(b)
reset(modelRea, data=myunemployment)
reset(modelRea, power=2, type="regressor", data=myunemployment)

## Harvey-Collier
harvtest(modelRea, order.by = ~ log(m/p), data=myunemployment)
harvtest(modelRea, order.by = ~ log(G), data=myunemployment)
harvtest(modelRea, order.by = ~ log(x), data=myunemployment)
harvtest(modelRea, data=myunemployment)

## Rainbow
raintest(modelRea, order.by = "mahalanobis", data=myunemployment)
}
\keyword{datasets}
