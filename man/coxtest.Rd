\name{coxtest}
\alias{coxtest}
\title{Cox Test for Comparing Non-Nested Models}

\description{
 \code{coxtest} performs the Cox test for comparing two non-nested models.
}

\usage{
coxtest(formula1, formula2, data = list())
}

\arguments{
  \item{formula1}{either a symbolic description for the first model to be tested,
    or a fitted object of class \code{"lm"}.}
  \item{formula2}{either a symbolic description for the second model to be tested,
    or a fitted object of class \code{"lm"}.}
  \item{data}{an optional data frame containing the variables in the
    model. By default the variables are taken from the environment
    which \code{coxtest} is called from.}
}

\details{
  The idea of the Cox test is the following: if the first model contains
  the correct set of regressors, then a fit of the regressors from the 
  second model to the fitted values from first model should have no further
  explanatory value. But if it has, it can be concluded that model 1 does not
  contain the correct set of regressors.
  
  Hence, to compare both models the fitted values of model
  1 are regressed on model 2 and vice versa. A Cox test statistic is computed
  for each auxiliary model which is asymptotically standard normally distributed.
  
  For further details, see the references.  
}

\value{
 An object of class \code{"anova"} which contains the estimate plus corresponding
 standard error, z test statistic and p value for each auxiliary test.
}

\references{
R. Davidson & J. MacKinnon (1981). Several Tests for Model Specification in the
Presence of Alternative Hypotheses. \emph{Econometrica}, \bold{49}, 781-793.

W. H. Greene (1993), \emph{Econometric Analysis}, 2nd ed. Macmillan Publishing Company, New York.

W. H. Greene (2003). \emph{Econometric Analysis}, 5th ed. New Jersey, Prentice Hall.
}

\seealso{\code{\link{jtest}}, \code{\link{encomptest}}}

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
