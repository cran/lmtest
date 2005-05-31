\name{grangertest}
\alias{grangertest}
\alias{grangertest.default}
\alias{grangertest.formula}
\title{Test for Granger Causality}
\description{
 \code{grangertest} is a generic function for performing 
 a test for Granger causality.
}
\usage{
\method{grangertest}{default}(x, y, order = 1, na.action = na.omit, \dots)
\method{grangertest}{formula}(formula, data = list(), \dots)
}
\arguments{
  \item{x}{either a bivariate series (in which case \code{y} has
    to be missing) or a univariate series of observations.}
  \item{y}{a univariate series of observations (if \code{x} is
    univariate, too).}
  \item{order}{integer specifying th order of lags to include
    in the auxiliary regression.}
  \item{na.action}{a function for eliminating \code{NA}s
    after aligning the series \code{x} and \code{y}.}
  \item{\dots}{further arguments passed to \code{\link{waldtest}}.}
  \item{formula}{a formula specification of a bivariate series like \code{y ~ x}.}
  \item{data}{an optional data frame containing the variables in the
    model. By default the variables are taken from the environment
    which \code{grangertest} is called from.}  
  }
\details{
Currently, the methods for the generic function \code{grangertest} only
perform tests for Granger causality in bivariate series. The test is simply
a Wald test comparing the unrestricted model---in which \code{y} is explained
by the lags (up to order \code{order}) of \code{y} and \code{x}---and the
restricted model---in which \code{y} is only explained by the lags of \code{y}.

Both methods are simply convenience interfaces to \code{\link{waldtest}}.
}

\value{
 An object of class \code{"anova"} which contains the residual degrees of freedom,
 the difference in degrees of freedom, Wald statistic and corresponding p value.
}

\seealso{\code{\link{waldtest}}, \code{\link{ChickEgg}}}

\examples{
## Which came first: the chicken or the egg?
data(ChickEgg)
grangertest(egg ~ chicken, order = 3, data = ChickEgg)
grangertest(chicken ~ egg, order = 3, data = ChickEgg)

## alternative ways of specifying the same test
grangertest(ChickEgg, order = 3)
grangertest(ChickEgg[, 1], ChickEgg[, 2], order = 3)
}

\keyword{htest}
