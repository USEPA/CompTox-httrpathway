library(openxlsx)
library(tcplfit2)
#' Pathway Concentration Response (P-value)
#'
#' Performs signature concentration response using p-value based cutoffs.
#'
#' Null dataset and dataset should have already been scored using signatureScore
#' and the given sigset and method. This function prepares signatureScore output
#' for CR processing, calls signatureConcRespCore_pval, formats the output,
#' saves it to disk, then calls plotouter for CR plots, if desired. If
#' conthits = T and nametag = NULL, the nametag "conthits" is automatically added
#' to the output file.
#'
#' @param sigset Name of the signature set.
#' @param dataset Name of the data set.
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#' @param nullset Name of the null data set.
#' @param conthits conthits = T uses continuous hitcalls, otherwise it's discrete.
#' @param aicc aicc = T uses corrected AIC to choose winning method; otherwise
#'   regular AIC.
#' @param minsigsize Minimum allowed signature size. Sample/signature combinations
#'   with less than this number of non-missing l2fc's will be discarded.
#'
#' @import data.table
#' @import parallel
#' @import openxlsx
#'
#' @return If to.file = T, nothing. If to.file = F, dataframe with signature CR
#'   output.
#' @export
signatureScoreCutoffDebug <- function(sigset="pilot_tiny",
                              sigcatalog="signatureDB_master_catalog 2020-05-05",
                              dataset="DMEM_6hr_pilot_normal_pe_1",
                              method="mygsea",
                              nullset="DMEM_6hr_pilot_normal_pe_1_RAND1000",
                              pval = .05,
                              conthits = T,
                              aicc = F,
                              minsigsize = 10) {

  printCurrentFunction(paste(dataset,sigset,method))
  starttime = proc.time()

  annotations <- signatureCatalogLoader(sigset,sigcatalog)

  #automatically add conthits nametag, then add _ to nametag if it exists
  nametag = "conthits"

  # get signaturescoremat
  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,".RData")
  print(file)
  load(file)
  #remove scores under path limit
  signaturescoremat = signaturescoremat[signaturescoremat$size >= minsigsize,]

  #attach pval cutoff and bmed
  pvalkey = getpvalcutoff(sigset, nullset, method, pvals = pval, numsds = 1)
  cat("> signatureConcResp 1\n")
  signaturescoremat$cutoff = pvalkey$cutoff[!is.na(pvalkey$pvalue)][match(signaturescoremat$signature, pvalkey$signature)]
  signaturescoremat$bmed = pvalkey$bmed[!is.na(pvalkey$pvalue)][match(signaturescoremat$signature, pvalkey$signature)]
  signaturescoremat$onesd = pvalkey$cutoff[!is.na(pvalkey$numsd)][match(signaturescoremat$signature, pvalkey$signature)]
  cat("> signatureConcResp 2\n")

  sig <- "CMAP fulvestrant 1e-06 100 8242 100"
  name <- "Fulvestrant"
  temp <- signaturescoremat[is.element(signaturescoremat$name,name),]
  nsig <- length(unique(signaturescoremat$signature))
  temp <- temp[is.element(temp$signature,sig),]
  cat(sigset,name,sig,nsig,temp[1,"cutoff"],temp[1,"bmed"],"\n")
  browser()
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
  browser()
}

