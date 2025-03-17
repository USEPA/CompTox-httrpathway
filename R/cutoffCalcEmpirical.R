#####################################################################################################
#' Calculate the signature-wise cutoffs based on the empirical distributions
#' which does not break any correlations between genes
#'
#' @param pval The p-value for the baseline distribution
#' @param nlowconc Only include the lowest nlowconc concentrations for each chemical
#' @param mc.cores Number of cores to use when running parallel
#' @param dtxsid.exclude dtxsids to exclude, default NULL
#' @param signaturescoremat Dataframe with one row per chemical/conc/signature combination. Columns
#'   are: sample_id, dtxsid, casrn, name, time, conc, sigset, signature, size, signature_score
#' @param CHEM_DICT Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#' @importFrom parallel makePSOCKcluster clusterExport parLapply stopCluster
#'
#' @return signature-wise cutoffs dataframe
#' @export cutoffCalcEmpirical
#####################################################################################################
cutoffCalcEmpirical <- function(
                      pval=0.05,
                      nlowconc=2,
                      mc.cores = 1,
                      dtxsid.exclude=NULL,
                      signaturescoremat,
                      CHEM_DICT){
  printCurrentFunction()

  chem_dict_filter <- CHEM_DICT[CHEM_DICT$conc_index<=nlowconc,]
  if(!is.null(dtxsid.exclude)) chem_dict_filter <- chem_dict_filter[!is.element(chem_dict_filter$dtxsid,dtxsid.exclude),]

  cdindex <- paste(chem_dict_filter$sample_id,chem_dict_filter$conc)
  signaturescoremat$index <- paste(signaturescoremat$sample_id,signaturescoremat$conc)
  SIGSCOREMAT <- signaturescoremat[is.element(signaturescoremat$index,cdindex),]

  siglist <- unique(SIGSCOREMAT$signature)
  cat("start calculating cutoffs\n")
  if(mc.cores > 1){
    cl <- makePSOCKcluster(mc.cores)
    #clusterExport(cl,c("SIGSCOREMAT","rowMean"))
    res <- parLapply(cl = cl, X=siglist, fun=cutoffCalc.inner.empirical, pval=pval, SIGSCOREMAT=SIGSCOREMAT)
    cutoffs <- do.call(rbind.data.frame,res)
  }
  else {
    res <- lapply(X=siglist, FUN=cutoffCalc.inner.empirical, pval, SIGSCOREMAT=SIGSCOREMAT)
    cutoffs <- do.call(rbind.data.frame,res)
  }
  if(mc.cores > 1) stopCluster(cl)
  cat("finish calculating cutoffs\n")
  return(cutoffs)
}
#####################################################################################################
#' Inner function for the cutoff calculation based on the empirical distributions
#'
#' @param signature The name of the signature for which the cutoff is to be calculated
#' @param pval The p-value for the baseline distribution
#' @param SIGSCOREMAT The signature score matrix
#' @importFrom stats qnorm median sd
#'
#' @return vector containing the signature, cutoff, sd (one standard deviation of the signature score matrix), bmed (median of the signature score matrix)
#' @export cutoffCalc.inner.empirical
#####################################################################################################
cutoffCalc.inner.empirical <- function(signature, pval, SIGSCOREMAT) {
  name.list <- c("signature","sd","bmed","cutoff","pvalue","numsd")
  cutoffs <- as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(cutoffs) <- name.list
  temp <- SIGSCOREMAT[is.element(SIGSCOREMAT$signature,signature),"signature_score"]

  cutoffs[1,"signature"] <- signature
  cutoffs[1,"sd"] <- 0
  cutoffs[1,"bmed"] <- 0
  cutoffs[1,"cutoff"] <- 0
  cutoffs[1,"pvalue"] <- pval
  cutoffs[1,"numsd"] <- 1

  scale <- 1-pval/2
  sdval <- sd(temp)
  cutoffs[1,"sd"] <- sdval
  cutoff <- qnorm(scale,mean=0,sd=sdval)
  cutoffs[1,"cutoff"] <- cutoff
  cutoffs[1,"bmed"] <- median(temp)
  return(as.vector(cutoffs))
}

