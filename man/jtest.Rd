\name{jtest}
\alias{jtest}
\title{J Test for Comparing Non-Nested Models}

\description{
 \code{jtest} performs the Davidson-MacKinnon J test for comparing non-nested
  models.
}

\usage{
 jtest(formula1, formula2, data = list(), vcov. = NULL, \dots)
}

\arguments{
  \item{formula1}{either a symbolic description for the first model to be tested,
    or a fitted object of class \code{"lm"}.}
  \item{formula2}{either a symbolic description for the second model to be tested,
    or a fitted object of class \code{"lm"}.}
  \item{data}{an optional data frame containing the variables in the
    model. By default the variables are taken from the environment
    which \code{jtest} is called from.}
  \item{vcov.}{a function for estimating the covariance matrix of the regression
    coefficients, e.g., \code{\link[sandwich]{vcovHC}}.}
  \item{\dots}{further arguments passed to \code{\link{coeftest}}.}
}

\details{
  The idea of the J test is the following: if the first model contains
  the correct set of regressors, then including the fitted values of the
  second model into the set of regressors should provide no significant 
  improvement. But if it does, it can be concluded that model 1 does not
  contain the correct set of regressors.
  
  Hence, to compare both models the fitted values of model
  1 are included into model 2 and vice versa. The J test statistic is simply
  the marginal test of the fitted values in the augmented model. This is
  performed by \code{\link{coeftest}}.
  
  For further details, see the references.  
}

\value{
 An object of class \code{"anova"} which contains the coefficient estimate
 of the fitted values in the augmented regression plus corresponding
 standard error, test statistic and p value.
}

\references{
R. Davidson & J. MacKinnon (1981). Several Tests for Model Specification in the
Presence of Alternative Hypotheses. \emph{Econometrica}, \bold{49}, 781-793.

W. H. Greene (1993), \emph{Econometric Analysis}, 2nd ed. Macmillan Publishing Company, New York.

W. H. Greene (2003). \emph{Econometric Analysis}, 5th ed. New Jersey, Prentice Hall.
}

\seealso{\code{\link{coxtest}}, \code{\link{encomptest}}}
\examples{
## Fit two competing, non-nested models for aggregate 
## consumption, as in Greene (1993), Examples 7.11 and 7.12

## load data and compute lags
data(USDistLag)
usdl <- na.contiguous(cbind(USDistLag, lag(USDistLag, k = -1)))
colnames(usdl) <- c("con", "gnp", "con1", "gnp1")

## C(t) = a0 + a1*Y(t) + a2*C(t-1) + u
fm1 <- lm(con ~ gnp + con1, data = usdl)

## C(t) = b0 + b1*Y(t) + b2*Y(t-1) + v
fm2 <- lm(con ~ gnp + gnp1, data = usdl)

## Cox test in both directions:
coxtest(fm1, fm2)

## ...and do the same for jtest() and encomptest().
## Notice that in this particular case they are coincident.
jtest(fm1, fm2)
encomptest(fm1, fm2)
}

\keyword{htest}
