\name{waldtest}
\alias{waldtest}
\alias{waldtest.formula}
\alias{waldtest.lm}
\title{Wald Test}

\description{
 \code{waldtest} is a generic function for carrying out Wald tests.
 Methods for comparing nested linear models (specified via \code{formula}s
 or fitted \code{"lm"} objects) are provided.
}

\usage{
 waldtest(object, \dots)

 \method{waldtest}{formula}(object, \dots)

 \method{waldtest}{lm}(object, \dots, vcov = NULL,
   test = c("F", "Chisq"), data = list()) 
}

\arguments{
  \item{object}{an object. For the \code{formula} and \code{lm} method
    this can either be a symbolic description of a linear model or
    a fitted object of class \code{"lm"}.}
  \item{\dots}{further arguments passed to methods. For the \code{formula}
    and \code{lm} method this can be further formulas or models to be
    compared with \code{object}. See below for further details.}
  \item{test}{character specifying wether to compute the finite sample
    F statistic (with approximate F distribution) or the large sample
    Chi-squared statistic (with asymptotic Chi-squared distribution).}
  \item{vcov}{a function for estimating the covariance matrix of the regression
    coefficients, e.g., \code{\link[sandwich]{vcovHC}}. If only two models
    are compared it can also be the covariance matrix of the more general
    model.}
  \item{data}{an optional data frame containing the variables in the
    model. By default the variables are taken from the environment
    which \code{waldtest} is called from.}
}

\details{
  \code{waldtest} is intended to be a generic function for comparisons
  of models via Wald tests. Methods for comparing linear models are provided.
  
  The methods for objects of classes \code{"formula"} and \code{"lm"} both
  compare linear models: \code{waldtest.formula} simply calls \code{waldtest.lm}
  passing on all other arguments. \code{waldtest.lm} takes all objects in
  \code{list(object, \dots)} and tries to compute fitted linear models from 
  them using the following algorithm: 
  \enumerate{
     \item If \code{object} is a formula, fit the corresponding linear model.
     \item For each two consecutive objects, \code{object1} and \code{object2}
           say, try to turn \code{object2} into a fitted model that can be
	   compared to \code{object1}.
     \item If \code{object2} is numeric, the corresponding element of
           \code{attr(terms(object1), "term.labels")} is selected.
     \item If \code{object2} is a character, the corresponding terms are
           included into an update formula like \code{. ~ . - term2a - term2b}.
     \item If \code{object2} is a formula, then compute the fitted model via
           \code{update(object1, object2)}.
  }
  Consequently, the models in \code{...} can be specified as integers, characters
  (both for terms that should be eliminated from the previous model), update formulas or
  fitted \code{"lm"} objects.
  
  Subsequently, a Wald test for each two consecutive models is carried out. This
  is very similar to \code{\link{anova}}, but with a few differences. If 
  only one \code{"lm"} object is specified, it is compared to the trivial model
  (with only an intercept). The test can be either the finite sample F statistic
  or the asymptotic Chi-squared statistic (\eqn{F = Chisq/k} if \eqn{k} is the difference
  in degrees of freedom). The covariance matrix is always estimated on the more general
  of two subsequent models (and not only in the most general model overall). If \code{vcov}
  is specified, HC and HAC estimators can also be plugged into \code{waldtest}.	     
}

\value{
 An object of class \code{"anova"} which contains the residual degrees of freedom
 in the encompassing model, the difference in degrees of freedom, Wald statistic
 (either \code{"F"} or \code{"Chisq"}) and corresponding p value.
}

\seealso{\code{\link{coeftest}}, \code{\link[stats]{anova}}, \code{\link[car]{linear.hypothesis}}}

\examples{
## fit two competing, non-nested models and their encompassing
## model for aggregate consumption, as in Greene (1993),
## Examples 7.11 and 7.12

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

## a simple ANOVA for fm3 vs. fm1
waldtest(fm3, fm2)
anova(fm3, fm2)
## as df = 1, the test is equivalent to the corresponding t test in
coeftest(fm3)

## various equivalent specifications of the two models
waldtest(fm3, fm2)
waldtest(fm3, 2)
waldtest(fm3, "con1")
waldtest(fm3, . ~ . - con1)

## comparing more than one model
## (euqivalent to the encompassing test)
waldtest(fm1, fm3, fm2)
encomptest(fm1, fm2)

## using the asymptotic Chisq statistic
waldtest(fm3, fm2, test = "Chisq")
## plugging in a HC estimator
if(require(sandwich)) waldtest(fm3, fm2, vcov = vcovHC)  
}

\keyword{htest}
