\name{hmctest}
\title{Harrison-McCabe-Test}
\usage{
hmctest(formula, T,  data=list())
}
\alias{hmctest}
\arguments{
 \item{formula}{a symbolic describtion for the model to be tested}
 \item{T}{split the model at time T}
 \item{data}{an optional data frame containing the variables in the model.
By default the variables are taken from the environment which 'hmctest' is
called from}
}
\description{
 \code{hmctest} performs the Harrison-McCabe-Test against heteroskedasticity.
}
\value{
 A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{not yet implemented, NA.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Kraemer, W., Sonnberger, H. (1986): The linear regression model
under test}
\author{Torsten Hothorn <Torsten.Hothorn@rzmail.uni-erlangen.de>}

\seealso{\code{\link{lm}}}

\examples{
x <- c(1:30);
err <- c(rnorm(10,0,1), rnorm(20,0,10));
y <- x^2 + err;
formular <- y ~ x;
hmc <- hmctest(formular, 10);
}
\keyword{linear models}