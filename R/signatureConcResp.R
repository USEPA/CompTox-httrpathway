#' Pathway Concentration Response (P-value)
#'
#' Performs signature concentration response using p-value based cutoffs.
#'
#' dataset should have already been scored using signatureScore and signatureScoreMerge
#' and the given sigset and method. This function prepares signatureScore output
#' for CR processing, calls tcplfit2::concRespCore, and formats the output.
#'
#' @param dataset Name of the data set.
#' @param sigset Name of the signature set.
#' @param sigcatalog Name of the signature catalog file
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param bmr_scale	bmr scaling factor. Default = 1.349
#' @param mc.cores Number of cores to parallelize with.
#' @param pval Desired cutoff p-value.
#' @param nlowconc Only include the lowest nlowconc concentrations for each chemical
#' @param aicc aicc = T uses corrected AIC to choose winning method; otherwise regular AIC
#' @param dtxsid.exclude dtxsids to exclude, default NULL
#' @param minsigsize Minimum allowed signature size. Sample/signature combinations
#'   with less than this number of non-missing l2fc's will be discarded.
#' @param fitmodels Vector of model names to use. Probably should include "cnst".
#' @param signaturescoremat dataframe returned by the upstream signatureScoreMerge function
#' @param sigdbgenelist full path to signature DB gene list file; default is repo version
#' @param CHEM_DICT Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#'
#' @importFrom data.table setDT
#' @importFrom parallel makePSOCKcluster clusterExport parLapplyLB stopCluster
#' @importFrom openxlsx read.xlsx
#' @import tcplfit2
#'
#' @return dataframe with signature CR output.
#' @export signatureConcResp
signatureConcResp <- function(dataset,
                              sigset,
                              sigcatalog="../inst/extdata/signatureDB_master_catalog_2022-05-16.xlsx",
                              method,
                              bmr_scale = 1.349,
                              mc.cores = 1,
                              pval = .05,
                              nlowconc = 2,
                              aicc = F,
                              dtxsid.exclude=NULL,
                              minsigsize = 10,
                              fitmodels = c("cnst", "hill", "poly1", "poly2", "pow", "exp2", "exp3", "exp4", "exp5"),
                              signaturescoremat,
                              sigdbgenelist="../inst/extdata/signatureDB_genelists.RDS",
                              CHEM_DICT) {

  printCurrentFunction(paste(dataset,sigset,method))
  starttime <- proc.time()

  annotations <- signatureCatalogLoader(sigset,sigcatalog,sigdbgenelist)

  cat("> signatureConcResp 0\n")
  #remove scores under path limit
  cat("> signatureConcResp 0.1\n")
  signaturescoremat <- signaturescoremat[signaturescoremat$size >= minsigsize,]

  # get the cutoffs
  pvalkey <- cutoffCalcEmpirical(signaturescoremat = signaturescoremat, CHEM_DICT = CHEM_DICT,
                                mc.cores = mc.cores, pval = pval, nlowconc = nlowconc, dtxsid.exclude = dtxsid.exclude)

  signaturescoremat <- as.data.table(signaturescoremat)

  cat("> signatureConcResp 1\n")
  signaturescoremat$cutoff <- pvalkey$cutoff[!is.na(pvalkey$pvalue)][match(signaturescoremat$signature, pvalkey$signature)]
  signaturescoremat$bmed <- pvalkey$bmed[!is.na(pvalkey$pvalue)][match(signaturescoremat$signature, pvalkey$signature)]
  signaturescoremat$onesd <- pvalkey$sd[!is.na(pvalkey$numsd)][match(signaturescoremat$signature, pvalkey$signature)]
  cat("> signatureConcResp 2\n")

  #cat("signaturescoremat is")
  #cat(signaturescoremat)

  #aggregate signaturescoremat by unique sample/signature per row; data table is considerably faster than aggregate
  if(is.element("time",names(signaturescoremat)))
    signaturescoremat <- signaturescoremat[, list(conc = list(conc), resp = list(signature_score), size = min(size)),
                                     by = c("sample_id", "dtxsid", "casrn", "name", "time", "signature", "bmed", "cutoff", "onesd")]
  else
    signaturescoremat <- signaturescoremat[, list(conc = list(conc), resp = list(signature_score), size = size),
                                                 by = c("sample_id", "dtxsid", "casrn", "name", "signature", "bmed", "cutoff", "onesd")]

  cat("> signatureConcResp 3\n")

  colnames(signaturescoremat)[colnames(signaturescoremat) == "size"] <- "signature_size"
  ordering <- order(tolower(signaturescoremat$name), tolower(signaturescoremat$signature))
  signaturescoremat <- signaturescoremat[ordering,]
  cat("> signatureConcResp 4\n")

  #turn signaturescoremat into a list of rows for lapply to use
  signaturescoremat <- as.list(as.data.frame(t(signaturescoremat), stringsAsFactors = F))
  cat("> signatureConcResp 5\n")

  if(mc.cores[1] > 1){
    cat("> signatureConcResp 6\n")
    cl <- makePSOCKcluster(mc.cores)
    clusterExport(cl, c("acy", "acgnlsobj", "bmdbounds", "bmdobj", "cnst", "exp2", "exp3", "exp4", "exp5", "fitcnst", "fithill", "fitgnls",
                        "fitcnst", "fitpoly1", "fitpoly2", "fitpow", "fitexp2", "fitexp3","fitexp4", "fitexp5",
                        "gnls" , "gnlsderivobj", "hillfn", "hitcontinner","hitloginner","tcplfit2_core","tcplhit2_core", "loggnls", "loghill", "nestselect",
                        "poly1", "poly2", "pow", "tcplObj", "toplikelihood"))

    SIGNATURE_CR <- parLapplyLB(cl=cl, X=signaturescoremat, fun=tcplfit2::concRespCore, fitmodels = fitmodels, bmr_scale=bmr_scale,
                             aicc=aicc, conthits=T, bmd_low_bnd=0.1, bmd_up_bnd=10, verbose=FALSE, chunk.size=ceiling(length(signaturescoremat)/5/mc.cores) )
    cat("> signatureConcResp 7\n")
  } else {
    cat("> signatureConcResp 6\n")
    SIGNATURE_CR <- lapply(X=signaturescoremat, FUN =tcplfit2::concRespCore, fitmodels = fitmodels, bmr_scale=bmr_scale,
                          aicc=aicc,conthits=T, bmd_low_bnd=0.1, bmd_up_bnd=10, verbose=F)
    cat("> signatureConcResp 7\n")
  }

  #construct SIGNATURE_CR
  SIGNATURE_CR <- as.data.frame(rbindlist(SIGNATURE_CR))
  rm(signaturescoremat)
  cat("> signatureConcResp 8\n")

  #add details to SIGNATURE_CR
  SIGNATURE_CR$target_class <- annotations$target_class[match(SIGNATURE_CR$signature, annotations$parent)]
  SIGNATURE_CR$super_target <- annotations$super_target[match(SIGNATURE_CR$signature, annotations$parent)]
  SIGNATURE_CR$sigset <- sigset
  SIGNATURE_CR$dataset <- dataset
  SIGNATURE_CR$method <- method
  SIGNATURE_CR$bmr_scale <- bmr_scale
  SIGNATURE_CR$ac50_loss <- as.numeric(SIGNATURE_CR$ac50_loss)
  cat("> signatureConcResp 9\n")


  cat("> signatureConcResp 10\n")

  if(mc.cores > 1) stopCluster(cl)
  print(proc.time() - starttime)
  return(SIGNATURE_CR)
}

