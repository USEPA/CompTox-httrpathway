#' Run All Pathway Concentration Response (P-Value)
#'
#' Driver for signature scoring and concentration response (CR).
#'
#' Signature scores are written to disk in output/signature_score_summary/.
#' Signature cutoffs are written to disk in output/signature_cutoff/.
#' CR results are written to disk in output/signature_conc_resp_summary/.
#'
#' @param dataset Name of data set.
#' @param sigset Name of signature set.
#' @param cutoff.dataset This is the data set name to sue when the cutoffs are taken from a different data set than
#'   the one currently being analyzed. The reason for doing this is if the current data set is small
#'   (small number of chemicals), and so not large enough to get a good estiamte of the underlying
#'   noise distribution. All of the other parameters for both data sets have to be the same
#' @param sigcatalog Name of the signature catalog
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param bmr_scale	bmr scaling factor. Default = 1.349
#' @param normfactor Factor to scale the native units up by to get onto a reasonable plotting value (~ -1 to 1)
#' @param minsigsize Minimum signature size.
#' @param pval P-value to use for noise estimation.
#' @param nlowconc Only include the lowest nlowconc concentrations for each chemical
#' @param mc.cores Vector with two values: number of cores to use for signature
#'   scoring and number of cores to use for CR. CR can usually handle the maximum
#'   number, but gsva scoring might require a smaller number to avoid memory
#'   overflow.
#' @param fitmodels Vector of model names to run conc/resp with. "cnst" should
#'   always be chosen.
#'
#' @return No output.
#'
#' remove gnls from default set
#' @export
runAllSignatureCR = function(dataset,
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
                             fitmodels = c("cnst", "hill",  "poly1", "poly2",
                                           "pow", "exp2", "exp3", "exp4", "exp5")){
  printCurrentFunction(paste(dataset,sigset,method))
  # load fcmat, chemdict
  basedir="../input/fcdata/"
  file <- paste0(basedir,"FCMAT2_",dataset,".RData")
  print(file)
  load(file)
  FCMAT2[is.nan(FCMAT2)] <- 0

  file <- paste0(basedir,"CHEM_DICT_",dataset,".RData",sep="")
  print(file)
  load(file)
  rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]

  cat("runAllSignatureCR: Start signatureScore\n")
  signatureScore(FCMAT2, CHEM_DICT, dataset=dataset, sigset=sigset,sigcatalog=sigcatalog,
                 method=method, normfactor=normfactor, mc.cores=mc.cores, minsigsize=minsigsize)

  cat("runAllSignatureCR: Start signatureScoreMerge\n")
  signatureScoreMerge(dataset=dataset, sigset=sigset, sigcatalog=sigcatalog, method=method)

  cat("runAllSignatureCR: Start cutoffCalcEmpirical\n")
  if(is.null(cutoff.dataset))
    cutoffCalcEmpirical(dataset=dataset, sigset=sigset, method=method, pval=pval,
                        nlowconc=nlowconc, mc.cores=mc.cores,do.load=T)

  cat("runAllSignatureCR: Start signatureConcResp\n")
  signatureConcResp(dataset=dataset, sigset=sigset, cutoff.dataset=cutoff.dataset, sigcatalog=sigcatalog,
                    method=method, pval=pval, nlowconc=nlowconc, minsigsize=minsigsize,
                    bmr_scale=bmr_scale, mc.cores=mc.cores, fitmodels=fitmodels)
}


