\name{custest}
\title{CUSUM-Test}
\usage{
raintest(formula, data=list())
}
\alias{custest}
\arguments{
 \item{formula}{a symbolic describtion for the model to be tested}
 \item{data}{an optional data frame containing the variables in the model.
By default the variables are taken from the environment which 'custest' is
called from}
}
\description{
 \code{custest} performs the CUSUM-Tests against structural change. This
test is based on recursive residuals and returns a random walk. 
}
\value{
 A vector containing the random-walk.
}
\references{Kraemer, W., Sonnberger, H. (1986): The linear regression model
under test}
\author{Torsten Hothorn <hothorn@amadeus.statistik.uni-dortmund.de>}

\seealso{\code{\link{lm}}}

\examples{
# generate regressor
x <- c(1:30);
# generate a structural change
y <- c(1*x[1:10], 3*x[11:30]) + rnorm(30,0,1);
# perform CUSUM-Test
rw <- custest(y ~x -1);
# plot random walk
plot(rw);
}
\keyword{linear models}