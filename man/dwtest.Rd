\name{dwtest}
\alias{dwtest}
\title{Durbin-Watson Test}
\description{
 Performs the Durbin-Watson test for autocorrelation of disturbances.
}
\usage{
dwtest(formula, iterations = 15, exact = NULL, tol = 1e-10,
  data = list())
}
\arguments{
 \item{formula}{a symbolic describtion for the model to be tested.}
 \item{iterations}{an integer specifying the number of iterations
    when calculating the p-value with the "pan" algorithm.}
 \item{exact}{logical. If set to \code{FALSE} a normal approximation
   will be used to compute the p value, if \code{TRUE} the "pan"
   algorithm is used. The default is to use "pan" if the sample size
   is < 100.}
 \item{tol}{tolerance. Eigenvalues computed have to be greater than
            \code{tol} to be treated as non-zero.}
 \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which \code{dwtest}
    is called from.}
}
\details{The Durbin-Watson test has the null hypothesis "rho = 0", this
implementation is built for the alternative "rho > 0".
The null distribution of the Durbin-Watson test statistic is a linear
combination of chi-squared distributions. The p value is computed using a
Fortran version of the Applied Statistics Algorithm AS 153 by Farebrother
(1980, 1984). This algorithm is called "pan" or "gradsol". For large sample
sizes the algorithm might fail to compute the p value; in that case a warning
is printed and an approximate p value will be given; this p value is computed
using a normal approximation with mean and variance of the Durbin-Watson test
statistic.

Examples can not only be found on this page, but also on the help pages of the
data sets \code{\link{bondyield}}, \code{\link{currencysubstitution}},
\code{\link{growthofmoney}}, \code{\link{moneydemand}}, \code{\link{unemployment}},
\code{\link{wages}}.
}
\value{An object of class \code{"htest"} containing:
\itemize
\item{statistic}{the test statistic.}
\item{p.value}{the corresponding p-value.}
\item{method}{a character string with the method used.}
\item{data.name}{a character string with the data name.}
}
\references{

J. Durbin & G.S. Watson (1950),
Testing for Serial Correlation in Least Squares Regression I.
\emph{Biometrika} \bold{37}, 409--428

J. Durbin & G.S. Watson (1951),
Testing for Serial Correlation in Least Squares Regression II.
\emph{Biometrika} \bold{38}, 159--178

J. Durbin & G.S. Watson (1971),
Testing for Serial Correlation in Least Squares Regression III.
\emph{Biometrika} \bold{58}, 1--19

R.W. Farebrother (1980),
Pan's Procedure for the Tail Probabilities of the
Durbin-Watson Statistic (Corr: 81V30 p189; AS R52: 84V33 p363- 366; AS
R53: 84V33 p366- 369).
\emph{Applied Statistics} \bold{29},  224--227

R. W. Farebrother (1984),
[AS R53] A Remark on Algorithms AS 106 (77V26 p92-98), AS 153 (80V29 p224-227)
  and AS 155: The Distribution of a Linear Combination of $\chi^2$ Random
  Variables (80V29 p323-333)
\emph{Applied Statistics} \bold{33}, 366--369

W. Kraemer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica


}

\seealso{\code{\link{lm}}}

\examples{

## generate two AR(1) error terms with parameter
## rho = 0 (white noise) and rho = 0.9 respectively
err1 <- rnorm(100)

## generate regressor and dependent variable
x <- rep(c(-1,1), 50)
y1 <- 1 + x + err1

## perform Durbin-Watson test
dwtest(y1 ~ x)

if(library(ts, logical = TRUE)) {
err2 <- filter(err1, 0.9, method="recursive")
y2 <- 1 + x + err2
dwtest(y2 ~ x) }

}
\keyword{htest}
