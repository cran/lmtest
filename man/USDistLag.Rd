\name{USDistLag}
\alias{USDistLag}
\title{US Macroeconomic Data}
\description{
US macroeconomic data for fitting a distributed lag model.
}
\usage{data(USDistLag)}
\format{
An annual time series from 1963 to 1982 with 2 variables.
 \describe{
    \item{consumption}{real consumption,}
    \item{gnp}{gross national product (deflated by CPI).}
 }
}

\source{Table 7.7 in Greene (1993)}

\references{
Greene W.H. (1993), \emph{Econometric Analysis}, 2nd edition.
  Macmillan Publishing Company, New York.

Executive Office of the President (1983), \emph{Economic Report of the
  President}. US Government Printing Office, Washington, DC.
}

\examples{
## Willam H. Greene, Econometric Analysis, 2nd Ed.
## Chapter 7
## load data set, p. 221, Table 7.7
data(USDistLag)

## fit distributed lag model, p.221, Example 7.8
usdl <- na.contiguous(cbind(USDistLag, lag(USDistLag, k = -1)))
colnames(usdl) <- c("con", "gnp", "con1", "gnp1")

## C(t) = b0 + b1*Y(t) + b2*C(t-1) + u
fm <- lm(con ~ gnp + con1, data = usdl)
summary(fm)
vcov(fm)
}

\keyword{datasets}
