#' Get P-Value Cutoff
#'
#' Retrieves pathway cutoffs for a given null dataset.
#'
#' Calculates median of all scores for a given pathway as well as a cutoff
#' based on the specified null dataset.
#' P-values represent the percentage of scores that are greater in distance
#' from the median than the cutoff. Numsd gives a cutoff that is the given
#' number of standard deviations from the median. Each row of the output
#' corresponds to one pathway and one pvalue or numsd. If both pvals and numsds
#' are specified, the output contains a column for each, and the unused
#' identifier(pvalue or numsd) in each row will contain NA.
#'
#' @param pathset Name of pathway set used to score null data.
#' @param nullset Name of null data set.
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#' @param pvals Vector of p-values to get cutoff for.
#' @param numsds Vector of number of standard deviations to get cutoff for.
#'   For instance, numsds = 1 will return cutoffs at 1 standard deviation.
#'
#' @import data.table
#'
#' @return Dataframe with 4 or 5 columns: pathway, cutoff, bmed (median of
#'   all samples for that pathway), pvalue (pvalue corresponding to each cutoff),
#'   numsd (number of sds corresponding to each cutoff).
#' @export
getpvalcutoff = function(pathset, nullset, method, pvals = NULL, numsds = NULL){
  printCurrentFunction()
  #get null pathway scores
  file <- paste0("../output/pathway_score_summary/PATHSCOREMAT_",pathset,"_",nullset,"_",method,".RData")
  load(file)

  #find p-value for each pathway
  if(!is.null(pvals)){
    pout = as.data.frame(setDT(pathscoremat)[, list(cutoff = quantile(abs(pathway_score-median(pathway_score)),1-pvals),
                                                   bmed = median(pathway_score)), by = list(pathway)])
    pout$pvalue = rep_len(pvals, nrow(pout))
    if(is.null(numsds)) output = pout else pout$numsd = NA_real_

  }

  #find x*standard deviation for each pathway
  if(!is.null(numsds)){
    sdout = as.data.frame(setDT(pathscoremat)[, list(cutoff = numsds*sd(pathway_score-median(pathway_score)) ,
                                                   bmed = median(pathway_score)), by = list(pathway)])
    sdout$numsd = rep_len(numsds, nrow(sdout))
    if(is.null(pvals)) output = sdout else {
      sdout$pvalue = NA_real_
      output = rbind(pout,sdout)
    }
  }

  return(output)
}
