\name{raintest}
\alias{raintest}
\encoding{latin1}
\title{Rainbow Test}

\description{Rainbow test for linearity.}

\usage{
raintest(formula, fraction = 0.5, order.by = NULL, center = NULL,
   data=list())
}

\arguments{
  \item{formula}{a symbolic description for the model to be tested
    (or a fitted \code{"lm"} object).}
 \item{fraction}{numeric. The percentage of observations in the subset is
   determined by
   \code{fraction*n} if \code{n} is the number of observations
   in the model.}
 \item{order.by}{Either a vector \code{z} or a formula with a single explanatory
    variable like \code{~ z}. The observations in the model
    are ordered by the size of \code{z}. If set to \code{NULL} (the
    default) the observations are assumed to be ordered (e.g., a
    time series). If set to \code{"mahalanobis"} then the observations
    are ordered by their Mahalanobis distances from the mean regressor.}
 \item{center}{numeric. If \code{center} is smaller than 1 it is
     interpreted as percentages of data, i.e. the subset is chosen that
     \code{n*fraction} observations are around observation number
     \code{n*center}. If \code{center} is greater than 1 it is interpreted to be
     the index of the center of the subset. By default \code{center} is 0.5.

     If the Mahalanobis distance is chosen \code{center} is taken to be
     the mean regressor, but can be specified to be a \code{k}-dimensional
     vector if \code{k} is the number of regressors and should be in the
     range of the respective regressors.}
\item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which
    \code{raintest} is called from.}
}

\details{The basic idea of the Rainbow test is that even if the true
 relationship is
non-linear, a good linear fit can be achieved on a subsample in the "middle" of
the data. The null hypothesis is rejected whenever the overall fit is
significantly worse than the fit for the subsample. The test statistic under
\eqn{H_0} follows an F distribution with \code{parameter} degrees of freedom.

Examples can not only be found on this page, but also on the help pages of the
data sets \code{\link{bondyield}}, \code{\link{currencysubstitution}},
\code{\link{growthofmoney}}, \code{\link{moneydemand}},
 \code{\link{unemployment}}, \code{\link{wages}}.
}

\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{parameter}{degrees of freedom.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{data.name}{a character string giving the name(s) of the data.}
}

\references{

J.M. Utts (1982),
The Rainbow Test for Lack of Fit in Regression.
\emph{Communications in Statistics -- Theory and Methods} \bold{11},
2801--2815.

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model under Test}. Heidelberg: Physica

}

\seealso{\code{\link{lm}}}

\examples{
x <- c(1:30)
y <- x^2 + rnorm(30,0,2)
rain <- raintest(y ~ x)
rain
## critical value
qf(0.95, rain$parameter[1], rain$parameter[2])
}

\keyword{htest}
