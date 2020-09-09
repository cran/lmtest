\name{coeftest}
\alias{coeftest}
\alias{coefci}
\alias{coeftest.default}
\alias{coeftest.survreg}
\alias{coeftest.glm}
\alias{coeftest.mlm}
\alias{coeftest.breakpointsfull}
\alias{print.coeftest}
\alias{confint.coeftest}
\alias{coef.coeftest}
\alias{df.residual.coeftest}
\alias{nobs.coeftest}
\alias{logLik.coeftest}
\alias{coefci.default}
\alias{coefci.survreg}
\alias{coefci.glm}
\alias{coefci.mlm}

\title{Inference for Estimated Coefficients}

\description{
  \code{coeftest} is a generic function for performing 
  z and (quasi-)t Wald tests of estimated coefficients.
  \code{coefci} computes the corresponding Wald confidence
  intervals. 
}

\usage{
coeftest(x, vcov. = NULL, df = NULL, \dots)

\method{coeftest}{default}(x, vcov. = NULL, df = NULL, \dots, save = FALSE)

coefci(x, parm = NULL, level = 0.95, vcov. = NULL, df = NULL, \dots)
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
  \item{\dots}{further arguments passed to the methods
    and to \code{vcov.} in the default method.}
  \item{save}{logical. Should the object \code{x} itself be
    saved as an attribute? (This may be useful for further
    processing of \code{coeftest} objects, e.g., as part of
    model summaries.)}
  \item{parm}{a specification of which parameters are to be given
    confidence intervals, either a vector of numbers or a vector
    of names.  If missing, all parameters are considered.}
  \item{level}{the confidence level required.}
}

\details{
The generic function \code{coeftest} currently has a default
method (which works in particular for \code{"lm"} objects) and
dedicated methods for objects of class
\code{"glm"} (as computed by \code{\link[stats]{glm}}),
\code{"mlm"} (as computed by \code{\link[stats]{lm}} with multivariate responses),
\code{"survreg"} (as computed by \code{\link[survival]{survreg}}), and
\code{"breakpointsfull"} (as computed by \code{\link[strucchange]{breakpoints.formula}}).

The default method assumes that a \code{coef} methods exists,
such that \code{coef(x)} yields the estimated coefficients.

To specify the corresponding covariance matrix \code{vcov.} to be used, there
are three possibilities:
1. It is pre-computed and supplied in argument \code{vcov.}.
2. A function for extracting the covariance matrix from 
\code{x} is supplied, e.g., \code{\link[sandwich]{sandwich}},
\code{\link[sandwich]{vcovHC}}, \code{\link[sandwich]{vcovCL}},
or \code{\link[sandwich]{vcovHAC}} from package \pkg{sandwich}.
3. \code{vcov.} is set to \code{NULL}, then it is assumed that
a \code{vcov} method exists, such that \code{vcov(x)} yields
a covariance matrix. Illustrations are provided in the examples below.

The degrees of freedom \code{df} determine whether a normal
approximation is used or a t distribution with \code{df} degrees
of freedom. The default method computes \code{df.residual(x)}
and if this is \code{NULL}, \code{0}, or \code{Inf} a z test is performed.
The method for \code{"glm"} objects always uses \code{df = Inf} (i.e., a z test).

The corresponding Wald confidence intervals can be computed either
by applying \code{coefci} to the original model or \code{\link[stats]{confint}}
to the output of \code{coeftest}. See below for examples.

Finally, \code{\link[stats]{nobs}} and \code{\link[stats]{logLik}}
methods are provided which work, provided that there are such methods
for the original object \code{x}. In that case, \code{"nobs"} and
\code{"logLik"} attributes are stored in the \code{coeftest} output
so that they can be still queried subsequently. If both methods are
available, \code{\link[stats]{AIC}} and \code{\link[stats]{BIC}}
can also be applied.
}

\value{
  \code{coeftest} returns an object of class \code{"coeftest"} which
  is essentially a coefficient matrix with columns containing the
  estimates, associated standard errors, test statistics and p values.
  Attributes for a \code{"method"} label, and the \code{"df"} are
  added along with \code{"nobs"} and \code{"logLik"} (provided that
  suitable extractor methods \code{\link[stats]{nobs}} and
  \code{\link[stats]{logLik}} are available). Optionally, the full
  object \code{x} can be \code{save}d in an attribute \code{"object"}
  to facilitate further model summaries based on the \code{coeftest}
  result.
  
  \code{coefci} returns a matrix (or vector) with columns giving
  lower and upper confidence limits for each parameter. These will
  be labeled as (1-level)/2 and 1 - (1-level)/2 in percent.
}

\seealso{\code{\link{lm}}, \code{\link{waldtest}}}

\examples{
## load data and fit model
data("Mandible", package = "lmtest")
fm <- lm(length ~ age, data = Mandible, subset=(age <= 28))

## the following commands lead to the same tests:
summary(fm)
(ct <- coeftest(fm))

## a z test (instead of a t test) can be performed by
coeftest(fm, df = Inf)

## corresponding confidence intervals
confint(ct)
coefci(fm)
## which in this simple case is equivalent to
confint(fm)

## extract further model information either from
## the original model or from the coeftest output
nobs(fm)
nobs(ct)
logLik(fm)
logLik(ct)
AIC(fm, ct)
BIC(fm, ct)

if(require("sandwich")) {
## a different covariance matrix can be also used:
(ct <- coeftest(fm, df = Inf, vcov = vcovHC))

## the corresponding confidence interval can be computed either as
confint(ct)
## or based on the original model
coefci(fm, df = Inf, vcov = vcovHC)

## note that the degrees of freedom _actually used_ can be extracted
df.residual(ct)
## which differ here from
df.residual(fm)

## vcov can also be supplied as a function with additional arguments
coeftest(fm, df = Inf, vcov = vcovHC, type = "HC0")
## or as a matrix
coeftest(fm, df = Inf, vcov = vcovHC(fm, type = "HC0"))
}
}

\keyword{htest}
