library(openxlsx)
library(tcplfit2)
#' Pathway Concentration Response (P-value)
#'
#' Performs signature concentration response using p-value based cutoffs.
#'
#' dataset should have already been scored using signatureScore
#' and the given sigset and method. This function prepares signatureScore output
#' for CR processing, calls signatureConcRespCore_pval, formats the output,
#' saves it to disk
#'
#' @param dataset Name of the data set.
#' @param sigset Name of the signature set.
#' @param cutoff.dataset This is the data set name to sue when the cutoffs are taken from a different data set than
#'   the one currently being analyzed. The reason for doing this is if the current data set is small
#'   (small number of chemicals), and so not large enough to get a good estiamte of the underlying
#'   noise distribution. All of the other parameters for both data sets have to be the same
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param bmr_scale	bmr scaling factor. Default = 1.349
#' @param mc.cores Number of cores to parallelize with.
#' @param pval Desired cutoff p-value.
#' @param nlowconc Only include the lowest nlowconc concentrations for each chemical
#' @param aicc aicc = T uses corrected AIC to choose winning method; otherwise
#'   regular AIC.
#' @param minsigsize Minimum allowed signature size. Sample/signature combinations
#'   with less than this number of non-missing l2fc's will be discarded.
#' @param fitmodels Vector of model names to use. Probably should include "cnst".
#'
#' @import data.table
#' @import parallel
#' @import openxlsx
#'
#' @return If to.file = T, nothing. If to.file = F, dataframe with signature CR
#'   output.
#' @export
signatureConcResp <- function(dataset,
                              sigset,
                              cutoff.dataset,
                              sigcatalog,
                              method,
                              bmr_scale = 1.349,
                              mc.cores = 1,
                              pval = .05,
                              nlowconc = 2,
                              aicc = F,
                              minsigsize = 10,
                              fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                            "exp4", "exp5")) {

  printCurrentFunction(paste(dataset,sigset,method))
  starttime = proc.time()

  annotations <- signatureCatalogLoader(sigset,sigcatalog)

  cat("> signatureConcResp 0\n")
  # get signaturescoremat
  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,".RData")
  load(file)
  #remove scores under path limit
  signaturescoremat = signaturescoremat[signaturescoremat$size >= minsigsize,]

  # get the cutoffs
  if(!is.null(cutoff.dataset))
    file = paste0("../output/signature_cutoff/signature_cutoff_",sigset,"_",cutoff.dataset,"_",method,"_",pval,"_",nlowconc,"_with_gene_correlations_empirical.xlsx")
  else
    file = paste0("../output/signature_cutoff/signature_cutoff_",sigset,"_",dataset,"_",method,"_",pval,"_",nlowconc,"_with_gene_correlations_empirical.xlsx")
  pvalkey = read.xlsx(file)

  cat("> signatureConcResp 1\n")
  signaturescoremat$cutoff = pvalkey$cutoff[!is.na(pvalkey$pvalue)][match(signaturescoremat$signature, pvalkey$signature)]
  signaturescoremat$bmed = pvalkey$bmed[!is.na(pvalkey$pvalue)][match(signaturescoremat$signature, pvalkey$signature)]
  signaturescoremat$onesd = pvalkey$sd[!is.na(pvalkey$numsd)][match(signaturescoremat$signature, pvalkey$signature)]
  cat("> signatureConcResp 2\n")

  #aggregate signaturescoremat by unique sample/signature per row; data table is considerably faster than aggregate
  if(is.element("time",names(signaturescoremat)))
    signaturescoremat = setDT(signaturescoremat)[, list(conc = list(conc),resp = list(signature_score), size = min(size)),
                                     by = list(sample_id, dtxsid, casrn, name, time, signature, bmed, cutoff, onesd)]
  else
    signaturescoremat = setDT(signaturescoremat)[, list(conc = list(conc),resp = list(signature_score), size = min(size)),
                                                 by = list(sample_id, dtxsid, casrn, name, signature, bmed, cutoff, onesd)]

  cat("> signatureConcResp 3\n")

  colnames(signaturescoremat)[colnames(signaturescoremat) == "size"] = "signature_size"
  ordering = order(tolower(signaturescoremat$name), tolower(signaturescoremat$signature))
  signaturescoremat = signaturescoremat[ordering,]
  cat("> signatureConcResp 4\n")

  #turn signaturescoremat into a list of rows for lapply to use
  signaturescoremat = as.list(as.data.frame(t(signaturescoremat), stringsAsFactors = F))
  cat("> signatureConcResp 5\n")

  if(mc.cores[1] > 1){
    cat("> signatureConcResp 6\n")
    cl = makePSOCKcluster(mc.cores)
    clusterExport(cl, c("acy", "acgnlsobj", "bmdbounds", "bmdobj", "cnst", "exp2", "exp3", "exp4", "exp5", "fitcnst", "fithill", "fitgnls",
                        "fitcnst", "fitpoly1", "fitpoly2", "fitpow", "fitexp2", "fitexp3","fitexp4", "fitexp5",
                        "gnls" , "gnlsderivobj", "hillfn", "hitcontinner","hitloginner","tcplfit2_core","tcplhit2_core", "loggnls", "loghill", "nestselect",
                        "poly1", "poly2", "pow", "tcplObj", "toplikelihood"))
    SIGNATURE_CR = parLapplyLB(cl=cl, X=signaturescoremat, fun=concRespCore, fitmodels = fitmodels, bmr_scale=bmr_scale,
                             aicc=aicc, conthits=T, bmd_low_bnd=0.1, bmd_up_bnd=10, verbose=FALSE, chunk.size=ceiling(length(signaturescoremat)/5/mc.cores) )
    cat("> signatureConcResp 7\n")
  } else {
    cat("> signatureConcResp 6\n")
    SIGNATURE_CR = lapply(X=signaturescoremat, FUN =concRespCore, fitmodels = fitmodels, bmr_scale=bmr_scale,
                          aicc=aicc,conthits=T, bmd_low_bnd=0.1, bmd_up_bnd=10, verbose=F)
    cat("> signatureConcResp 7\n")
  }

  #construct SIGNATURE_CR
  SIGNATURE_CR = as.data.frame(rbindlist(SIGNATURE_CR))
  rm(signaturescoremat)
  cat("> signatureConcResp 8\n")

  #add details to SIGNATURE_CR
  SIGNATURE_CR$target_class <- annotations$target_class[match(SIGNATURE_CR$signature, annotations$parent)]
  SIGNATURE_CR$super_target <- annotations$super_target[match(SIGNATURE_CR$signature, annotations$parent)]
  SIGNATURE_CR$sigset = sigset
  SIGNATURE_CR$dataset = dataset
  SIGNATURE_CR$method = method
  SIGNATURE_CR$bmr_scale = bmr_scale
  SIGNATURE_CR$ac50_loss = as.numeric(SIGNATURE_CR$ac50_loss)
  cat("> signatureConcResp 9\n")

  #save SIGNATURE_CR
  dir.create("../output/signature_conc_resp_summary/", showWarnings = F)
  file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_", pval, "_conthits.RData")
  if(bmr_scale!=1.349) file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_bmr_scale_",bmr_scale,"_", pval, "_conthits.RData")
  save(SIGNATURE_CR,file=file)

  cat("> signatureConcResp 10\n")

  if(mc.cores > 1) stopCluster(cl)
  print(proc.time() - starttime)
}

