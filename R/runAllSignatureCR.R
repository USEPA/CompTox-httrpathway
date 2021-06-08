#' Run All Pathway Concentration Response (P-Value)
#'
#' Driver for signature scoring and concentration response (CR).
#'
#' CR requires signature scores to have already been computed for a nullset.
#' randomdata() can generate a nullset, and this function can compute signature
#' scores for it by setting dataset = nullset and nullset = NULL. Pathway
#' scores are written to disk in output/signature_score_summary/. CR results
#' are written to disk in output/signature/conc_resp_summary/.
#'
#' @param basedir Folder that stores FCMAT2 and CHEM_DICT files.
#' @param dataset Name of data set.
#' @param sigset Name of signature set.
#' @param sigcatalog Name of the signature catalog
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param bmr_scale	bmr scaling factor. Default = 1.349
#' @param normfactor Factor to scale the native units up by to get onto a reasonable plotting value (~ -1 to 1)
#' @param minsigsize Minimum signature size.
#' @param conthits conthits = T uses continous hitcall; conthits = F uses discrete
#'   hitcalls.
#' @param nullset Name of null dataset. Set nullset = NULL to skip CR.
#' @param do.plot do.plot = T generates a CR plot for every sample/signature
#'   combination.
#' @param do.cr Run the concentration-response step (set to FALSE for the null set)
#' @param pval P-value to use for noise estimation.
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
runAllSignatureCR = function(basedir="../input/fcdata/",
                             dataset="DMEM_6hr_pilot_normal_pe_1",
                             sigset,
                             sigcatalog,
                             method = "gsea",
                             bmr_scale=1.349,
                             normfactor=7500,
                             minsigsize = 10,
                             conthits = T,
                             nullset,
                             do.plot = T,
                             do.cr=T,
                             pval = .05,
                             mc.cores = c(1,1),
                             fitmodels = c("cnst", "hill",  "poly1", "poly2",
                                           "pow", "exp2", "exp3", "exp4", "exp5")){
  printCurrentFunction(paste(dataset,sigset,method,nullset))
  # load fcmat, chemdict
  file <- paste0(basedir,"FCMAT2_",dataset,".RData")
  print(file)
  load(file)
  FCMAT2[is.nan(FCMAT2)] <- 0

  file <- paste0(basedir,"CHEM_DICT_",dataset,".RData",sep="")
  print(file)
  load(file)
  rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]

  #run signature scores
  cat("runAllSignatureCR: Start signatureScore\n")
  signatureScore(FCMAT2, CHEM_DICT,sigset,sigcatalog,dataset,method=method,
                 normfactor=normfactor,
                 mc.cores=mc.cores[1], minsigsize = minsigsize)
  if(is.null(nullset)) return()

  cat("runAllSignatureCR: Start signatureScoreMerge\n")
  signatureScoreMerge(sigset,sigcatalog,dataset,method,nullset)

  if(do.cr) {

    #run conc/resp
    cat("runAllSignatureCR: Start signatureConcResp\n")
    signatureConcResp(sigset,sigcatalog,dataset,method=method, bmr_scale=bmr_scale, nullset = nullset, mc.cores=mc.cores[2], do.plot = do.plot,
                      to.file=T, pval = pval, minsigsize = minsigsize, conthits = conthits,
                      fitmodels = fitmodels)
  }

}

