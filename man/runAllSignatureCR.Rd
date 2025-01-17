% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runAllSignatureCR.R
\name{runAllSignatureCR}
\alias{runAllSignatureCR}
\title{Run All Pathway Concentration Response (P-Value)}
\usage{
runAllSignatureCR(
  dataset,
  sigset,
  cutoff.dataset,
  sigcatalog,
  method,
  bmr_scale = 1.349,
  normfactor = 7500,
  minsigsize = 10,
  pval = 0.05,
  nlowconc = 2,
  mc.cores = 1,
  fitmodels = c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3", "exp4",
    "exp5")
)
}
\arguments{
\item{dataset}{Name of data set.}

\item{sigset}{Name of signature set.}

\item{cutoff.dataset}{This is the data set name to sue when the cutoffs are taken from a different data set than
the one currently being analyzed. The reason for doing this is if the current data set is small
(small number of chemicals), and so not large enough to get a good estiamte of the underlying
noise distribution. All of the other parameters for both data sets have to be the same}

\item{sigcatalog}{Name of the signature catalog}

\item{method}{Pathway scoring method in c("fc", "gsva", "gsea")}

\item{bmr_scale}{bmr scaling factor. Default = 1.349}

\item{normfactor}{Factor to scale the native units up by to get onto a reasonable plotting value (~ -1 to 1)}

\item{minsigsize}{Minimum signature size.}

\item{pval}{P-value to use for noise estimation.}

\item{nlowconc}{Only include the lowest nlowconc concentrations for each chemical}

\item{mc.cores}{Vector with two values: number of cores to use for signature
scoring and number of cores to use for CR. CR can usually handle the maximum
number, but gsva scoring might require a smaller number to avoid memory
overflow.}

\item{fitmodels}{Vector of model names to run conc/resp with. "cnst" should
always be chosen.}
}
\value{
No output.

remove gnls from default set
}
\description{
Driver for signature scoring and concentration response (CR).
}
\details{
Signature scores are written to disk in output/signature_score_summary/.
Signature cutoffs are written to disk in output/signature_cutoff/.
CR results are written to disk in output/signature_conc_resp_summary/.
}
