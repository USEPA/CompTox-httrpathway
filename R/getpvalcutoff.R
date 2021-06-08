#' Get P-Value Cutoff
#'
#' Retrieves signature cutoffs for a given null dataset.
#'
#' Calculates median of all scores for a given signature as well as a cutoff
#' based on the specified null dataset.
#' P-values represent the percentage of scores that are greater in distance
#' from the median than the cutoff. Numsd gives a cutoff that is the given
#' number of standard deviations from the median. Each row of the output
#' corresponds to one signature and one pvalue or numsd. If both pvals and numsds
#' are specified, the output contains a column for each, and the unused
#' identifier(pvalue or numsd) in each row will contain NA.
#'
#' @param pathset Name of signature set used to score null data.
#' @param nullset Name of null data set.
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param pvals Vector of p-values to get cutoff for.
#' @param numsds Vector of number of standard deviations to get cutoff for.
#'   For instance, numsds = 1 will return cutoffs at 1 standard deviation.
#' @param verbose If TRUE, write extra output
#' @import data.table
#'
#' @return Dataframe with 4 or 5 columns: signature, cutoff, bmed (median of
#'   all samples for that signature), pvalue (pvalue corresponding to each cutoff),
#'   numsd (number of sds corresponding to each cutoff).
#' @export
getpvalcutoff = function(pathset, nullset, method, pvals = NULL, numsds = NULL,verbose=T){
  printCurrentFunction(paste(pathset,":",nullset))

  #get null signature scores
  file <- paste0("../output/signature_score_summary/signaturescoremat_",pathset,"_",nullset,"_",method,".RData")
  if(verbose) print(file)
  load(file)
  if(verbose) cat("  file loaded\n")

  ###########################################
  ###########################################
  debug=F
  if(debug) {
    sig <- "CMAP fulvestrant 1e-06 100 8242 100"
    temp <- signaturescoremat[is.element(signaturescoremat$signature,sig),]
    mval <- mean(temp$signature_score)
    cat(pathset,nrow(temp),mean(temp$signature_score),sd(temp$signature_score),"\n")
    browser()
  }
  ###########################################
  ###########################################
  if(verbose) cat("  find p-value for each signature\n")
  if(!is.null(pvals)){
    pout = as.data.frame(setDT(signaturescoremat)[, list(cutoff = quantile(abs(signature_score-median(signature_score)),1-pvals,na.rm=T),
                                                   bmed = median(signature_score)), by = list(signature)])
    pout$pvalue = rep_len(pvals, nrow(pout))
    if(is.null(numsds)) output = pout else pout$numsd = NA_real_

  }

  if(verbose) cat("  find x*standard deviation for each signature\n")
  if(!is.null(numsds)){
    sdout = as.data.frame(setDT(signaturescoremat)[, list(cutoff = numsds*sd(signature_score-median(signature_score)) ,
                                                   bmed = median(signature_score)), by = list(signature)])
    sdout$numsd = rep_len(numsds, nrow(sdout))
    if(is.null(pvals)) output = sdout else {
      sdout$pvalue = NA_real_
      output = rbind(pout,sdout)
    }
  }

  return(output)
}
