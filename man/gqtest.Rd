\name{gqtest}
\title{Goldfeld-Quandt-Test}
\usage{
gqtest(formula, T, data=list())
}
\alias{gqtest}
\arguments{
 \item{formula}{a symbolic describtion for the model to be tested}
 \item{T}{split the model at time T}
 \item{data}{an optional data frame containing the variables in the model.
By default the variables are taken from the environment which 'gqtest' is
called from}
}
\description{
 \code{gqtest} performs the Goldfeld-Quandt-Test against heteroskadasticity.
Given a known time \code{T}, the Goldfeld-Quandt tests the null-hypothesis:
variances at time 1..T and T..n are equal.
}
\details{
The test statistic under H0 follows a F distribution with df1 and df2 degree
of freedom.
}
\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Kraemer, W., Sonnberger, H. (1986): The linear regression model
under test}
\author{Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de>}

\seealso{\code{\link{lm}}}

\examples{
x <- c(1:30)
err <- c(rnorm(10,0,1), rnorm(20,0,10))
y <- x + err
formular <- y ~ x
gq <- gqtest(formular, 10)
qf(0.95, gq$parameter[1], gq$parameter[2])
}
\keyword{htest}