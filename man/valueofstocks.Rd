\name{valueofstocks}
\alias{valueofstocks}
\docType{data}
\encoding{latin1}
\title{ Value of Stocks  }
\usage{data(valueofstocks)}
\description{
Value of Stocks Data
}
\format{
 A multivariate quarterly time series from 1960(1) to 1977(3) with variables
 \describe{
    \item{VST}{value of stocks,}
    \item{MB}{monetary base,}
    \item{RTPD}{dollar rent on producer durables,}
    \item{RTPS}{dollar rent on producer structures,}
    \item{XBC}{production capacity for business output.}
  }

}
\source{The data was originally studied by Woglom (1981), the data set is given
 in
Krämer and Sonnberger (1986).
}
\references{

G. Woglom (1981),
A Reexamination of the Role of Stocks in the Consumption Function and the
Transmission Mechanism.
\emph{Journal of Money, Credit and Banking} \bold{13}, 215--220

W. Krämer & H. Sonnberger (1986),
\emph{The Linear Regression Model Under Test}. Heidelberg: Physica


}
\examples{
data(valueofstocks)
lm(log(VST) ~., data=valueofstocks)
}
\keyword{datasets}
