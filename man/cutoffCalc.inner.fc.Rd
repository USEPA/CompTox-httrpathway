% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cutoffCalc.R
\name{cutoffCalc.inner.fc}
\alias{cutoffCalc.inner.fc}
\title{Inner function for the cutoff calculation based on the analytical method}
\usage{
cutoffCalc.inner.fc(parent, catalog, allgenes, pval)
}
\arguments{
\item{parent}{The name of the signature parent for which the cutoff is to be calculated}

\item{catalog}{The signature catalog}

\item{allgenes}{The list of all the genes in the data set}

\item{pval}{The p-value for the baseline distribution}
}
\value{
vector containing the parent (signature), cutoff, sd, bmed
}
\description{
Inner function for the cutoff calculation based on the analytical method
}
