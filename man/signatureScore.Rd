% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signatureScore.R
\name{signatureScore}
\alias{signatureScore}
\title{Signature Score}
\usage{
signatureScore(
  FCMAT2,
  CHEM_DICT,
  sigset,
  sigcatalog = "../inst/extdata/signatureDB_master_catalog_2022-05-16.xlsx",
  method,
  normfactor = 7500,
  mc.cores = 1,
  minsigsize = 10,
  sigdbgenelist = "../inst/extdata/signatureDB_genelists.RDS"
)
}
\arguments{
\item{FCMAT2}{Sample by gene matrix of log2(fold change)'s. Rownames are
sample keys and colnames are genes.}

\item{CHEM_DICT}{Dataframe with one row per sample key and seven columns:
sample_key, sample_id, conc, time, casrn, name, dtxsid.}

\item{sigset}{Name of signature set.}

\item{sigcatalog}{full path to signature catalog xlsx file; default is repo version}

\item{method}{Signature scoring method in c("fc", "gsva", "gsea")}

\item{normfactor}{Value passed to the plotting code to scale the y values}

\item{mc.cores}{Number of cores to use.}

\item{minsigsize}{Minimum allowed signature size BEFORE accounting for
missing values.}

\item{sigdbgenelist}{full path to signature DB gene list file; default is repo version}
}
\value{
Returns data frame of signature scores
}
\description{
Computes and saves signature scores.
}
\details{
signatureScore is a driver for various scoring methods. The three that are
currently available are "gsva", "gsea", "fc".
Beware running out of memory on large runs with gsva, Linux, and many cores -- ensure your system has enough memory allocated depending on data size.
Signature size is counted according to number of genes in the signature that are
also in the column names of FCMAT2. However, each method performs a more
rigorous size count internally that accounts for missing values and adds this
to the output.
}
