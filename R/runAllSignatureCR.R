#' Run All Pathway Concentration Response (P-Value)
#'
#' Driver for signature scoring and concentration response (CR).
#'
#'
#' @param dataset Name of data set.
#' @param sigset Name of signature set.
#' @param sigcatalog Name of the signature catalog
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param bmr_scale benchmark response (bmr) scaling factor. Default = 1.349
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
#' @param FCMAT2 Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param CHEM_DICT Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#' @param sigdbgenelist full path to signature DB gene list file; default is repo version
#'
#' @return dataframe with signature concentration-response results
#'
#' @export runAllSignatureCR
runAllSignatureCR = function(dataset,
                             sigset,
                             sigcatalog="../inst/extdata/signatureDB_master_catalog_2022-05-16.xlsx",
                             method,
                             bmr_scale = 1.349,
                             normfactor = 7500,
                             minsigsize = 10,
                             pval = 0.05,
                             nlowconc = 2,
                             mc.cores = 1,
                             fitmodels = c("cnst", "hill",  "poly1", "poly2", "pow", "exp2", "exp3", "exp4", "exp5"),
                             FCMAT2,
                             CHEM_DICT,
                             sigdbgenelist="../inst/extdata/signatureDB_genelists.RDS"){
  printCurrentFunction(paste(dataset,sigset,method))
  # load fcmat, chemdict
  FCMAT2[is.nan(FCMAT2)] <- 0

  rownames(CHEM_DICT) <- CHEM_DICT$sample_key

  cat("runAllSignatureCR: Start signatureScore\n")
  signaturescoremat <- signatureScore(FCMAT2, CHEM_DICT, sigset=sigset,sigcatalog=sigcatalog,
                 method=method, normfactor=normfactor, mc.cores=mc.cores, minsigsize=minsigsize, sigdbgenelist=sigdbgenelist)

  cat("runAllSignatureCR: Start signatureScoreMerge\n")
  signaturescoremat <- signatureScoreMerge(sigset=sigset, sigcatalog=sigcatalog, method=method, signaturescoremat=signaturescoremat, sigdbgenelist=sigdbgenelist)

  # cat("runAllSignatureCR: Start cutoffCalcEmpirical\n")
  # cutoffCalcEmpirical(pval=pval, nlowconc=nlowconc, mc.cores=mc.cores, signaturescoremat=signaturescoremat, CHEM_DICT=CHEM_DICT) #this is already nested within the signatureConcResp function

  cat("runAllSignatureCR: Start signatureConcResp\n")
  return(signatureConcResp(dataset=dataset, sigset=sigset,sigcatalog=sigcatalog,
                    method=method, pval=pval, nlowconc=nlowconc, minsigsize=minsigsize,
                    bmr_scale=bmr_scale, mc.cores=mc.cores, fitmodels=fitmodels,signaturescoremat=signaturescoremat, sigdbgenelist=sigdbgenelist, CHEM_DICT=CHEM_DICT))
}


