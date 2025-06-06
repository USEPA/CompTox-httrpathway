% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statfuns.R
\name{WRMSE}
\alias{WRMSE}
\title{Weighted Root-mean-square-error}
\usage{
WRMSE(x, y, w)
}
\arguments{
\item{x}{First vector of numbers.}

\item{y}{Second vector of numbers.}

\item{w}{Vector of weights.}
}
\value{
Weighted RMSE.
}
\description{
Computes root-mean-square error with weighted average.
}
\details{
x,y,w should all be the same length. Order of x and y won't change output.
}
\examples{
WRMSE(1:3, c(1,3,5), 1:3)
}
