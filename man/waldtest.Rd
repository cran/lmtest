\name{waldtest}
\alias{waldtest}
\alias{waldtest.formula}
\alias{waldtest.default}
\alias{waldtest.lm}
\title{Wald Test of Nested Models}

\description{
 \code{waldtest} is a generic function for carrying out Wald tests.
 The default method can be employed for comparing nested (generalized)
 linear models (see details below).
}

\usage{
 waldtest(object, \dots)

 \method{waldtest}{default}(object, \dots, vcov = NULL,
   test = c("Chisq", "F"), name = NULL) 

 \method{waldtest}{formula}(object, \dots, data = list())

 \method{waldtest}{lm}(object, \dots, test = c("F", "Chisq"))
}

\arguments{
  \item{object}{an object. See below for details.}
  \item{\dots}{further object specifications passed to methods. See below for details.}
  \item{vcov}{a function for estimating the covariance matrix of the regression
    coefficients, e.g., \code{\link[sandwich]{vcovHC}}. If only two models
    are compared it can also be the covariance matrix of the more general
    model.}
  \item{test}{character specifying whether to compute the large sample
    Chi-squared statistic (with asymptotic Chi-squared distribution) or the finite sample
    F statistic (with approximate F distribution).}
  \item{name}{a function for extracting a suitable name/description from
    a fitted model object. By default the name is queried by calling
    \code{\link{formula}}.}
  \item{data}{a data frame containing the variables in the model.}
}

\details{
  \code{waldtest} is intended to be a generic function for comparisons
  of models via Wald tests. The default method consecutively compares
  the fitted model object \code{object} with the models passed in \code{\dots}.
  Instead of passing the fitted model objects in \code{\dots}, several other
  specifications are possible. For all objects in \code{list(object, \dots)}
  the function tries to consecutively compute fitted models using the following
  updating algorithm: 
  \enumerate{
     \item For each two consecutive objects, \code{object1} and \code{object2}
           say, try to turn \code{object2} into a fitted model that can be
	   compared to (the already fitted model object) \code{object1}.
     \item If \code{object2} is numeric, the corresponding element of
           \code{attr(terms(object1), "term.labels")} is selected to be omitted.
     \item If \code{object2} is a character, the corresponding terms are
           included into an update formula like \code{. ~ . - term2a - term2b}.
     \item If \code{object2} is a formula, then compute the fitted model via
           \code{update(object1, object2)}.
  }
  Consequently, the models in \code{\dots} can be specified as integers, characters
  (both for terms that should be eliminated from the previous model), update formulas or
  fitted model objects. Except for the last case, the existence of an \code{\link{update}}
  method is assumed. See also the examples for an illustration.
  
  Subsequently, a Wald test for each two consecutive models is carried out. This
  is similar to \code{\link{anova}} (which typically performs likelihood-ratio tests),
  but with a few differences. If only one fitted model object is specified, it is compared
  to the trivial model (with only an intercept). The test can be either the finite sample
  F statistic or the asymptotic Chi-squared statistic (\eqn{F = Chisq/k} if \eqn{k} is the
  difference in degrees of freedom). The covariance matrix is always estimated on the more general
  of two subsequent models (and not only in the most general model overall). If \code{vcov}
  is specified, HC and HAC estimators can also be plugged into \code{waldtest}.	     

  The default method is already very general and applicable to a broad range of fitted
  model objects, including \code{\link{lm}} and \code{\link{glm}} objects. It can be
  easily made applicable to other model classes as well by providing suitable methods
  to the standard generics \code{\link{terms}} (for determining the variables in the model
  along with their names), \code{\link{update}} (unless only fitted model objects are passed
  to \code{waldtest}, as mentioned above), \code{\link{nobs}} (or \code{\link{residuals}}, used for determining
  the number of observations), \code{\link{df.residual}} (needed only for the F statistic),
  \code{\link{coef}} (for extracting the coefficients; needs to be named matching the names
  in \code{terms}), \code{\link{vcov}} (can be user-supplied; needs to be named matching the
  names in \code{terms}). Furthermore, some means of determining a suitable \code{name} for
  a fitted model object can be specified (by default this is taken to be the result of
  a call to \code{\link{formula}}, if available).
  
  The \code{"formula"} method fits a \code{\link{lm}} first and then calls the \code{"lm"}
  method. The \code{"lm"} method just calls the default method, but sets the default
  test to be the F test.
}

\value{
 An object of class \code{"anova"} which contains the residual degrees of freedom,
 the difference in degrees of freedom, Wald statistic
 (either \code{"Chisq"} or \code{"F"}) and corresponding p value.
}

\seealso{\code{\link{coeftest}}, \code{\link[stats]{anova}}, \code{\link[car]{linearHypothesis}}}

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

## a simple ANOVA for fm3 vs. fm2
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
## (equivalent to the encompassing test)
waldtest(fm1, fm3, fm2)
encomptest(fm1, fm2)

## using the asymptotic Chisq statistic
waldtest(fm3, fm2, test = "Chisq")
## plugging in a HC estimator
if(require(sandwich)) waldtest(fm3, fm2, vcov = vcovHC)  
}

\keyword{htest}
