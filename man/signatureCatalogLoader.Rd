% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signatureCatalogLoader.R
\name{signatureCatalogLoader}
\alias{signatureCatalogLoader}
\title{Merge the up and down halves of the pathway data}
\usage{
signatureCatalogLoader(
  sigset = "wgcna",
  sigcatalog = "signatureDB_wgcna_mcf7_ph1_pe1_normal_good_pg_MCF7_12_10_catalog"
)
}
\arguments{
\item{sigset}{Name of the signature set.}

\item{sigcatlog}{Nmae of the catalog file}
}
\value{
the trimmed signature table
}
\description{
Merge the up and down halves of the pathway data
}
