% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cutoffCalcEmpirical.R
\name{cutoffCalcEmpirical}
\alias{cutoffCalcEmpirical}
\title{Calculate the signature-wise cutoffs based on the empirical distributions
which does not break any correlations between genes}
\usage{
cutoffCalcEmpirical(
  basedir = "../input/fcdata/",
  dataset = "heparg2d_toxcast_pfas_pe1_normal",
  sigset = "screen_large",
  method = "fc",
  pval = 0.05,
  nlowconc = 2,
  mc.cores = 1,
  dtxsid.exclude = NULL,
  do.load = T
)
}
\arguments{
\item{basedir}{Directory that holds FCMAT2 and CHEM_DICT files.}

\item{dataset}{Name of actual dataset to base cutoff on.}

\item{sigset}{THe signature set}

\item{method}{The scoring method, either fc or gsea}

\item{pval}{The p-value for the baseline distribution}

\item{nlowconc}{Only include the lowest nlowconc concentrations for each chemical}

\item{mc.cores}{NUmber of cores to use when running parallel}

\item{dtxsid.exclude}{dtxsids to exclude, default NULL}

\item{do.load}{If TRUE, reload the FCMAT2 matrix, signature catalog and chemical dictionary, and store in globals}

\item{sigcatalog}{The name of the signature catalog to use}

\item{do.cov}{If TRUE, calculate the covariance matrix and store in a global}

\item{do.compare}{If TRUE, compare the cutoffs with those from the original method with no gene-gene correlation}

\item{to.file}{If TRUE, and do.compare=TRUE, send a plot of the comparison to a file}

\item{verbose}{If TRUE, write a line for each signature to show progress.}
}
\value{
No output.
}
\description{
Calculate the signature-wise cutoffs based on the empirical distributions
which does not break any correlations between genes
}
