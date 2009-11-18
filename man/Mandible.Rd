\name{Mandible}
\alias{Mandible}
\title{Mandible Data}
\usage{data(Mandible)}
\description{
Mandible Data.
}
\format{
  Data from 167 fetuses, especially:
  \describe{
    \item{age}{gestational age in weeks.}
    \item{length}{mandible length in mm.}
  }
}
\source{
  The data was originally published by Chitty et al., 1993, and analyzed in
Royston and Altman, 1994 (the data is given there). Only measurements with
\code{age <= 28} were used in this analysis.
}
\references{

  L. S. Chitty and S. Campbell and D. G. Altman (1993), 
  Measurement of the fetal mandible -- feasibility and
  construction of a centile chart., \emph{Prenatal Diagnosis}, 
  \bold{13}, 749--756.

  P. Royston and D. G. Altman (1994), Regression Using 
  Fractional Polynomials of Continuous Covariates: Parsimonious 
  Parametric Modelling. \emph{Applied Statistics}, \bold{43}, 429--453.

}

\examples{
data(Mandible)
lm(length ~ age, data=Mandible, subset=(age <= 28))

}
\keyword{datasets}
