\name{jocci}
\alias{fyff}
\alias{gmdc}
\alias{ip}
\alias{jocci}
\alias{lhur}
\alias{pw561}
\title{U.S. Macroeconomic Time Series}
\usage{
data(fyff)
data(gmdc)
data(ip)
data(jocci)
data(lhur)
data(pw561)
}
\description{
Several macroeconomic time series from the U.S.
}
\format{
 All data sets are multivariate monthly time series from
 1959(8) to 1993(12) (except 1993(10) for \code{jocci}) with variables
 \describe{
    \item{y}{original time series,}
    \item{dy}{transformed times series (first differences
              or log first differences),}
    \item{dy1}{transformed series at lag 1,}
    \item{dy2}{transformed series at lag 2,}
    \item{dy3}{transformed series at lag 3,}
    \item{dy4}{transformed series at lag 4,}
    \item{dy5}{transformed series at lag 5,}
    \item{dy6}{transformed series at lag 6.}
  }
}

\details{
The description from Stock & Watson (1996) for the time series
(with the transformation used):
  \describe{
    \item{fyff}{interest rate (first differences),}
    \item{gmdc}{pce, implicit price deflator: pce (1987 = 100) (log first differences),}
    \item{ip}{index of industrial production (log first differences),}
    \item{jocci}{department of commerce commodity price index (log first differences),}
    \item{lhur}{unemployment rate: all workers, 16 years & over (\%, sa) (first differences),}
    \item{pw561}{producer price index: crude petroleum (82 = 100, nsa) (log first differences).}
  }
Stock & Watson (1996) fitted an AR(6) model to all transformed time series.
}

\source{
Stock & Watson (1996) study the stability of 76 macroeconomic time series, which can be obtained
from Mark W. Watson's homepage at \url{http://www.princeton.edu/~mwatson/ddisk/bivtvp.zip}.
}
\references{

J.H. Stock & M.W. Watson (1996),
Evidence on Structural Instability in Macroeconomic Time Series Relations.
\emph{Journal of Business & Economic Statistics} \bold{14}, 11--30.

}
\examples{
data(jocci)

dwtest(dy ~ 1, data = jocci)
bgtest(dy ~ 1, data = jocci)
ar6.model <- dy ~ dy1 + dy2 + dy3 + dy4 + dy5 +dy6
bgtest(ar6.model, data = jocci)

var.model <- ~ I(dy1^2) + I(dy2^2) + I(dy3^2) + I(dy4^2) + I(dy5^2) + I(dy6^2)
bptest(ar6.model, var.model, data = jocci)
}
\keyword{datasets}
