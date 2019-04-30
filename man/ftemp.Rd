\name{ftemp}
\alias{ftemp}
\title{Femal Temperature Data}
\usage{data(ftemp)}
\description{
Daily morning temperature of adult female (in degrees Celsius).
}
\format{
  Univariate daily time series of 60 observations starting from 1990-07-11.
}
\details{
The data gives the daily morning temperature of an adult woman
measured in degrees Celsius at about 6.30am each morning.

At the start of the period the woman was sick, hence the high temperature.
Then the usual monthly cycle can be seen. On the second cycle, the temperature
doesn't complete the downward part of the pattern due to a conception.
}

\source{
The data set was taken from the Time Series Data Library, maintained by Rob Hyndman.
}

\examples{
data(ftemp)
plot(ftemp)
y <- window(ftemp, start = 8, end = 60)
if(require(strucchange)) {
  bp <- breakpoints(y ~ 1)
  plot(bp)
  fm.seg <- lm(y ~ 0 + breakfactor(bp))
  plot(y)
  lines(8:60, fitted(fm.seg), col = 4)
  lines(confint(bp))
}

}
\keyword{datasets}
