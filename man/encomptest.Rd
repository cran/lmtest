\name{encomptest}
\alias{encomptest}
\title{Encompassing Test for Comparing Non-Nested Models}

\description{
 \code{encomptest} performs the encompassing test of Davidson & MacKinnon
 for comparing non-nested models.
}

\usage{
 encomptest(formula1, formula2, data = list(), vcov. = NULL, \dots)
}

\arguments{
  \item{formula1}{either a symbolic description for the first model to be tested,
    or a fitted object of class \code{"lm"}.}
  \item{formula2}{either a symbolic description for the second model to be tested,
    or a fitted object of class \code{"lm"}.}
  \item{data}{an optional data frame containing the variables in the
    model. By default the variables are taken from the environment
    which \code{encomptest} is called from.}
  \item{vcov.}{a function for estimating the covariance matrix of the regression
    coefficients, e.g., \code{\link[sandwich]{vcovHC}}.}
  \item{\dots}{further arguments passed to \code{\link{waldtest}}.}
}

\details{
  To compare two non-nested models, the encompassing test fits an encompassing
  model which contains all regressors from both models such that the two
  models are nested within the encompassing model. A Wald test for comparing
  each of the models with the encompassing model is carried out by \code{\link{waldtest}}.
  
  For further details, see the references.  
}

\value{
 An object of class \code{"anova"} which contains the residual degrees of freedom
 in the encompassing model, the difference in degrees of freedom, Wald statistic
 (either \code{"F"} or \code{"Chisq"}) and corresponding p value.
}

\references{
R. Davidson & J. MacKinnon (1993). \emph{Estimation and Inference in Econometrics}.
New York, Oxford University Press.

W. H. Greene (1993), \emph{Econometric Analysis}, 2nd ed. Macmillan Publishing Company, New York.

W. H. Greene (2003). \emph{Econometric Analysis}, 5th ed. New Jersey, Prentice Hall.
}

\seealso{\code{\link{coxtest}}, \code{\link{jtest}}}
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

## Encompassing model
fm3 <- lm(con ~ gnp + con1 + gnp1, data = usdl)

## Cox test in both directions:
coxtest(fm1, fm2)

## ...and do the same for jtest() and encomptest().
## Notice that in this particular case they are coincident.
jtest(fm1, fm2)
encomptest(fm1, fm2)

## the encompassing test is essentially
waldtest(fm1, fm3, fm2)
}

\keyword{htest}
