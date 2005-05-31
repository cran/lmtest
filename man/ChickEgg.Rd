\name{ChickEgg}
\alias{ChickEgg}
\title{Chickens, Eggs, and Causality}
\description{
US chicken population and egg production.
}
\usage{data(ChickEgg)}
\format{
An annual time series from 1930 to 1983 with 2 variables.
 \describe{
    \item{chicken}{number of chickens (December 1 population
      of all US chickens excluding commercial broilers),}
    \item{egg}{number of eggs (US egg production in millions
      of dozens).}
 }
}

\source{The data set was provided by Walter Thurman and made available
for R by Roger Koenker. Unfortunately, the data is slightly different
than the data analyzed in Thurman & Fisher (1988).}

\references{
Thurman W.N. & Fisher M.E. (1988), Chickens, Eggs, and Causality, or
Which Came First?, \emph{American Journal of Agricultural Economics},
237-238.
}

\examples{
## Which came first: the chicken or the egg?
data(ChickEgg)
## chickens granger-cause eggs?
grangertest(egg ~ chicken, order = 3, data = ChickEgg)
## eggs granger-cause chickens?
grangertest(chicken ~ egg, order = 3, data = ChickEgg)

## To perform the same tests `by hand', you can use dynlm() and waldtest():
if(require(dynlm)) {
  ## chickens granger-cause eggs?
  em <- dynlm(egg ~ L(egg, 1) + L(egg, 2) + L(egg, 3), data = ChickEgg)
  em2 <- update(em, . ~ . + L(chicken, 1) + L(chicken, 2) + L(chicken, 3))
  waldtest(em, em2)

  ## eggs granger-cause chickens?
  cm <- dynlm(chicken ~ L(chicken, 1) + L(chicken, 2) + L(chicken, 3), data = ChickEgg)
  cm2 <- update(cm, . ~ . + L(egg, 1) + L(egg, 2) + L(egg, 3))
  waldtest(cm, cm2)
}
}

\keyword{datasets}
