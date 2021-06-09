#####################################################################################################
#' Calculate the signature-wise cutoffs based on the empirical distributions
#' which does not break any correlations between genes
#'
#' @param basedir Directory that holds FCMAT2 and CHEM_DICT files.
#' @param dataset Name of actual dataset to base cutoff on.
#' @param sigcatalog The name of the signature catalog to use
#' @param sigset THe signature set
#' @param method The scoring method, either fc or gsea
#' @param pval The p-value for the baseline distribution
#' @param nlowconc Only include the lowest nlowconc concentrations for each chemical
#' @param mc.cores NUmber of cores to use when running parallel
#' @param dtxsid.exclude dtxsids to exclude, default NULL
#' @param do.load If TRUE, reload the FCMAT2 matrix, signature catalog and chemical dictionary, and store in globals
#' @param do.cov If TRUE, calculate the covariance matrix and store in a global
#' @param do.compare If TRUE, compare the cutoffs with those from the original method with no gene-gene correlation
#' @param to.file If TRUE, and do.compare=TRUE, send a plot of the comparison to a file
#' @param verbose If TRUE, write a line for each signature to show progress.
#'
#' @return No output.
#' @export
#####################################################################################################
cutoffCalcEmpirical = function(basedir="../input/fcdata/",
                      dataset="heparg2d_toxcast_pfas_pe1_normal",
                      sigset="screen_large",
                      method="fc",
                      pval=0.05,
                      nlowconc=2,
                      mc.cores = 1,
                      dtxsid.exclude=NULL,
                      do.load=T){
  printCurrentFunction()

  if(do.load) {
    #load chem_dict
    file = paste0(basedir,"CHEM_DICT_",dataset,".RData")
    load(file)
    chem_dict_filter = CHEM_DICT[CHEM_DICT$conc_index<=nlowconc,]
    if(!is.null(dtxsid.exclude)) chem_dict_filter = chem_dict_filter[!is.element(chem_dict_filter$dtxsid,dtxsid.exclude),]

    cdindex = paste(chem_dict_filter$sample_id,chem_dict_filter$conc)
    #load signaturescoremat_
    file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,".RData")
    cat("  file:",file,"\n")
    load(file=file)
    signaturescoremat$index = paste(signaturescoremat$sample_id,signaturescoremat$conc)
    signaturescoremat$index = paste(signaturescoremat$sample_id,signaturescoremat$conc)
    SIGSCOREMAT <<- signaturescoremat[is.element(signaturescoremat$index,cdindex),]
  }
  siglist = unique(SIGSCOREMAT$signature)
  cat("start calculating cutoffs\n")
  if(mc.cores > 1){
    cl = makePSOCKcluster(mc.cores)
    clusterExport(cl,c("SIGSCOREMAT","rowMean"))
    res = parLapply(cl = cl, X=siglist, fun=cutoffCalc.inner.emperical, pval=pval)
    cutoffs = do.call(rbind,res)
  }
  else {
    res = lapply(X=siglist,FUN=cutoffCalc.inner.emperical,pval)
    cutoffs = do.call(rbind,res)
  }
  if(mc.cores > 1) stopCluster(cl)
  cat("finish calculating cutoffs\n")
  file = paste0("../output/signature_cutoff/signature_cutoff_",sigset,"_",dataset,"_",method,"_",pval,"_",nlowconc,"_with_gene_correlations_empirical.xlsx")
  write.xlsx(cutoffs,file)
}
#####################################################################################################
#' Inner function for the cutoff calculation
#'
#' @param signature The name of the signature  for which the cutoff is to be calculated
#' @param pval The p-value for the baseline distribution
#' @param covmat THe covariance matrix
#'
#' @return vector containing the signature, cutoff, sd, bmed
#' @export
#####################################################################################################
cutoffCalc.inner.emperical = function(signature, pval) {
  name.list = c("signature","sd","bmed","cutoff","pvalue","numsd")
  cutoffs = as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(cutoffs) = name.list
  temp = SIGSCOREMAT[is.element(SIGSCOREMAT$signature,signature),"signature_score"]

  cutoffs[1,"signature"] = signature
  cutoffs[1,"sd"] = 0
  cutoffs[1,"bmed"] = 0
  cutoffs[1,"cutoff"] = 0
  cutoffs[1,"pvalue"] = pval
  cutoffs[1,"numsd"] = 1

  scale = 1-pval/2
  sdval = sd(temp)
  cutoffs[1,"sd"] = sdval
  cutoff = qnorm(scale,mean=0,sd=sdval)
  cutoffs[1,"cutoff"] = cutoff
  cutoffs[1,"bmed"] = median(temp)
  return(as.vector(cutoffs))
}

