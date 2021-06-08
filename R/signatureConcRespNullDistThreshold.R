#' Export the Null set thresholds
#'
#' Null dataset and dataset should have already been scored using signatureScore
#' and the given sigset and method.
#' @param sigset Name of the signature set.
#' @param dataset Name of the data set.
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param pval Desired cutoff p-value.
#'
signatureConcRespNullDistThreshold <- function(sigset="screen_large",
                                               dataset="mcf7_ph1_pe1_normal_all_pg",
                                               method="fc",
                                               pval = .05) {
  printCurrentFunction(paste(dataset,sigset,method))
  nullset = paste0(dataset,"_RAND1000")

  cat("> start getting pval\n")
  pvalkey = getpvalcutoff(sigset, nullset, method, pvals = pval, numsds = 1,verbose=T)
  res1 = pvalkey[!is.na(pvalkey$pvalue),1:4]
  res1$type="pval"
  res2 = pvalkey[!is.na(pvalkey$numsd),c(1,2,3,5)]
  res2$type="numsd"
  names(res1)[4] = "value"
  names(res2)[4] = "value"
  res = rbind(res1,res2)
  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,"_nullSetThreshold.RData")
  signature_null_thresholds = res
  save(signature_null_thresholds,file=file)
}

