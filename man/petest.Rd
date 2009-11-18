\name{petest}
\alias{petest}
\title{PE Test for Linear vs. Log-Linear Specifications}

\description{
 \code{petest} performs the MacKinnon-White-Davidson PE test for comparing
  linear vs. log-linear specifications in linear regressions.
}

\usage{
 petest(formula1, formula2, data = list(), vcov. = NULL, \dots)
}

\arguments{
  \item{formula1}{either a symbolic description for the first model to be tested,
    or a fitted object of class \code{"lm"}.}
  \item{formula2}{either a symbolic description for the second model to be tested,
    or a fitted object of class \code{"lm"}.}
  \item{data}{an optional data frame containing the variables in the
    model. By default the variables are taken from the environment
    which \code{petest} is called from.}
  \item{vcov.}{a function for estimating the covariance matrix of the regression
    coefficients, e.g., \code{\link[sandwich]{vcovHC}}.}
  \item{\dots}{further arguments passed to \code{\link{coeftest}}.}
}

\details{
  The PE test compares two non-nest models where one has a linear
  specification of type \code{y ~ x1 + x2} and the other has a log-linear
  specification of type \code{log(y) ~ z1 + z2}. Typically, the
  regressors in the latter model are logs of the regressors in the
  former, i.e., \code{z1} is \code{log(x1)} etc.
  
  The idea of the PE test is the following: If the linear specification is
  correct then adding an auxiliary regressor with the difference of
  the log-fitted values from both models should be non-significant.
  Conversely, if the log-linear specification is correct then adding
  an auxiliary regressor with the difference of fitted values in levels
  should be non-significant. The PE test statistic is simply the marginal
  test of the auxiliary variable(s) in the augmented model(s). In \code{petest}
  this is performed by \code{\link{coeftest}}.
  
  For further details, see the references.  
}

\value{
 An object of class \code{"anova"} which contains the coefficient estimate
 of the auxiliary variables in the augmented regression plus corresponding
 standard error, test statistic and p value.
}

\references{
W.H. Greene (2003). \emph{Econometric Analysis}, 5th edition.
Upper Saddle River, NJ: Prentice Hall.

J. MacKinnon, H. White, R. Davidson (1983). Tests for Model Specification in the
Presence of Alternative Hypotheses: Some Further Results.
\emph{Journal of Econometrics}, \bold{21}, 53-70.

M. Verbeek (2004). \emph{A Guide to Modern Econometrics}, 2nd ed. Chichester, UK: John Wiley.
}

\seealso{\code{\link{jtest}}, \code{\link{coxtest}}, \code{\link{encomptest}}}
\examples{
if(require("AER")) {
## Verbeek (2004), Section 3
data("HousePrices", package = "AER")

### Verbeek (2004), Table 3.3
hp_lin <- lm(price ~ . , data = HousePrices)
summary(hp_lin)

### Verbeek (2004), Table 3.2
hp_log <- update(hp_lin, log(price) ~ . - lotsize + log(lotsize))
summary(hp_log)

## PE test
petest(hp_lin, hp_log)


## Greene (2003), Example 9.8
data("USMacroG", package = "AER")

## Greene (2003), Table 9.2
usm_lin <- lm(m1 ~ tbill + gdp, data = USMacroG)
usm_log <- lm(log(m1) ~ log(tbill) + log(gdp), data = USMacroG)
petest(usm_lin, usm_log)
## matches results from Greene's errata
}
}

\keyword{htest}
