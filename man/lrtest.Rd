\name{lrtest}
\alias{lrtest}
\alias{lrtest.formula}
\alias{lrtest.default}
\title{Likelihood Ratio Test of Nested Models}

\description{
 \code{lrtest} is a generic function for carrying out likelihood ratio tests.
 The default method can be employed for comparing nested (generalized)
 linear models (see details below).
}

\usage{
 lrtest(object, \dots)

 \method{lrtest}{default}(object, \dots, name = NULL) 

 \method{lrtest}{formula}(object, \dots, data = list())
}

\arguments{
  \item{object}{an object. See below for details.}
  \item{\dots}{further object specifications passed to methods. See below for details.}
  \item{name}{a function for extracting a suitable name/description from
    a fitted model object. By default the name is queried by calling
    \code{\link{formula}}.}
  \item{data}{a data frame containing the variables in the model.}
}

\details{
  \code{lrtest} is intended to be a generic function for comparisons
  of models via asymptotic likelihood ratio tests. The default method consecutively compares
  the fitted model object \code{object} with the models passed in \code{\dots}.
  Instead of passing the fitted model objects in \code{\dots}, several other
  specifications are possible. The updating mechanism is the same as for \code{\link{waldtest}}:
  the models in \code{\dots} can be specified as integers, characters
  (both for terms that should be eliminated from the previous model), update formulas or
  fitted model objects. Except for the last case, the existence of an \code{\link{update}}
  method is assumed. See \code{\link{waldtest}} for details.
  
  Subsequently, an asymptotic likelihood ratio test for each two consecutive models is carried out:
  Twice the difference in log-likelihoods (as derived by the \code{\link{logLik}} methods)
  is compared with a Chi-squared distribution.

  The \code{"formula"} method fits a \code{\link{lm}} first and then calls the default 
  method.
}

\value{
 An object of class \code{"anova"} which contains the log-likelihood, degrees of freedom,
 the difference in degrees of freedom, likelihood ratio Chi-squared statistic and corresponding p value.
}

\seealso{\code{\link{waldtest}}}

\examples{
## with data from Greene (1993):
## load data and compute lags
data("USDistLag")
usdl <- na.contiguous(cbind(USDistLag, lag(USDistLag, k = -1)))
colnames(usdl) <- c("con", "gnp", "con1", "gnp1")

fm1 <- lm(con ~ gnp + gnp1, data = usdl)
fm2 <- lm(con ~ gnp + con1 + gnp1, data = usdl)

## various equivalent specifications of the LR test
lrtest(fm2, fm1)
lrtest(fm2, 2)
lrtest(fm2, "con1")
lrtest(fm2, . ~ . - con1)
}

\keyword{htest}
