\name{coeftest}
\alias{coeftest}
\alias{coeftest.default}
\alias{coeftest.survreg}
\alias{coeftest.glm}
\alias{coeftest.mlm}
\alias{coeftest.breakpointsfull}
\alias{print.coeftest}
\title{Testing Estimated Coefficients}
\description{
 \code{coeftest} is a generic function for performing 
   z and (quasi-)t tests of estimated coefficients.
}
\usage{
coeftest(x, vcov. = NULL, df = NULL, \dots)
}
\arguments{
  \item{x}{an object (for details see below).}
  \item{vcov.}{a specification of the covariance
    matrix of the estimated coefficients. This can be
    specified as a matrix or as a function yielding
    a matrix when applied to \code{x}.}
  \item{df}{the degrees of freedom to be used. If this
    is a finite positive number a t test with \code{df}
    degrees of freedom is performed. In all other cases,
    a z test (using a normal approximation) is performed.
    By default it tries to use \code{x$df.residual}
    and performs a z test if this is \code{NULL}.}
  \item{\dots}{further arguments passed to the methods.}
  }
\details{
The generic function \code{coeftest} currently has a default
method (which works in particular for \code{"lm"} and
\code{"glm"} objects) and a method for objects of class
\code{"breakpointsfull"}
(as computed by \code{\link[strucchange]{breakpoints.formula}}).

The default method assumes that a \code{coef} methods exists,
such that \code{coef(x)} yields the estimated coefficients.

To specify a covariance matrix \code{vcov.} to be used, there
are three possibilities:
1. It is pre-computed and supplied in argument \code{vcov.}.
2. A function for extracting the covariance matrix from 
\code{x} is supplied, e.g., \code{\link[sandwich]{vcovHC}}
or \code{\link[sandwich]{vcovHAC}} from package \pkg{sandwich}.
3. \code{vcov.} is set to \code{NULL}, then it is assumed that
a \code{vcov} method exists, such that \code{vcov(x)} yields
a covariance matrix. For illustrations see below.

The degrees of freedom \code{df} determine whether a normal
approximation is used or a t distribution with \code{df} degrees
of freedoms is used. The default method uses \code{df.residual(x)}
and if this is \code{NULL} a z test is performed.
}

\value{
  An object of class \code{"coeftest"} which is essentially
  a coefficient matrix with columns containing the estimates,
  associated standard errors, test statistics and p values.
}

\seealso{\code{\link{lm}}, \code{\link{waldtest}}}

\examples{
## load data and fit model
data(Mandible)
fm <- lm(length ~ age, data=Mandible, subset=(age <= 28))

## the following commands lead to the same tests:
summary(fm)
coeftest(fm)

## a z test (instead of a t test) can be performed by
coeftest(fm, df = Inf)

if(require(sandwich)) {
## a different covariance matrix can be also used:
## either supplied as a function
coeftest(fm, df = Inf, vcov = vcovHC)
## or as a matrix
coeftest(fm, df = Inf, vcov = vcovHC(fm, type = "HC0"))
}
}

\keyword{htest}
