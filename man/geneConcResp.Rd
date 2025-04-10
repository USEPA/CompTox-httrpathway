% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneConcResp.R
\name{geneConcResp}
\alias{geneConcResp}
\title{Gene Concentration Response}
\usage{
geneConcResp(
  mc.cores = 20,
  to.file.path = NULL,
  pval = 0.05,
  aicc = F,
  fitmodels = c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3", "exp4", "exp5"),
  genefile = NULL,
  FCMAT2,
  CHEM_DICT
)
}
\arguments{
\item{mc.cores}{Number of parallel cores to use.}

\item{to.file.path}{when provided, path of RDS file where results are written to}

\item{pval}{P-value cutoff between 0 and 1.}

\item{aicc}{If aicc = T, corrected AIC is used instead of first order
(regular) AIC.}

\item{fitmodels}{Vector of models names to be used. Default is all of them.}

\item{genefile}{An optional file (.xlsx) that can be used to filter concentration-response modeling for a subset of genes of interest}

\item{FCMAT2}{chem/conc by gene or chem/conc by probe. Uses two lowest concentration of each
column to estimate noise cutoff (as opposed to signature CR).}

\item{CHEM_DICT}{Dataframe with one row per sample key and seven columns:
sample_key, sample_id, conc, time, casrn, name, dtxsid.}
}
\value{
dataframe of concentration response modeling results
}
\description{
Wrapper that performs concentration response modeling for gene/probe log2 fold-change values
}
\details{
Uses two lowest concentration of each
column to estimate noise cutoff (as opposed to signature CR).
}
