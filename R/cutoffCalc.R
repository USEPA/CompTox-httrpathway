#####################################################################################################
#' Calculate the signature-wise cutoffs based on the analytical method
#' which does not break any correlations between genes
#' and compare with already defined signature_cr object (see exportSignatureCutoffs function)
#'
#' @param sigcatalog The name of the signature catalog to use
#' @param sigset The signature set
#' @param pval The p-value for the baseline distribution
#' @param seed Random seed
#' @param nlowconc Only include the nlowconc concentrations for each chemical when calcualting the cutoff
#' @param mc.cores Number of cores to use when running parallel
#' @param dtxsid.exclude dtxsids to exclude, default NULL
#' @param do.cov If TRUE, calculate the covariance matrix and store in a global
#' @param do.compare If TRUE, compare the cutoffs with those from the original method with no gene-gene correlation
#' @param to.file If TRUE, and do.compare=TRUE, send a plot of the comparison to a file
#' @param sigdbgenelist full path to the signatureDB_genelist.RDS file
#' @param FCMAT2 Sample by gene matrix of log2(fold change)'s. Rownames are
#'   sample keys and colnames are genes.
#' @param CHEM_DICT Dataframe with one row per sample key and seven columns:
#'   sample_key, sample_id, conc, time, casrn, name, dtxsid.
#' @importFrom parallel makePSOCKcluster clusterExport parLapply stopCluster
#' @importFrom graphics par plot lines
#' @importFrom grDevices pdf dev.off
#' @importFrom openxlsx write.xlsx
#' @importFrom stats cov
#'
#' @return calculated cutoffs dataframe
#' @export cutoffCalc
#####################################################################################################
cutoffCalc = function(sigcatalog,
                      sigset,
                      pval=0.05,
                      seed = 12345,
                      nlowconc=2,
                      mc.cores = 1,
                      dtxsid.exclude=NULL,
                      do.cov=T,
                      do.compare=F,
                      to.file=F,
                      sigdbgenelist,
                      FCMAT2,
                      CHEM_DICT){
  printCurrentFunction()
  set.seed(seed)

  chem_dict_filter = CHEM_DICT[CHEM_DICT$conc_index<=nlowconc,]
  if(!is.null(dtxsid.exclude)) chem_dict_filter = chem_dict_filter[!is.element(chem_dict_filter$dtxsid,dtxsid.exclude),]

  sk.list = chem_dict_filter$sample_key
  FCMAT2[is.nan(FCMAT2)] = 0
  FCMAT2[is.na(FCMAT2)] = 0
  FCMAT2_lowconc <<- FCMAT2[sk.list,]

  file = paste0(sigcatalog)
  catalog = read.xlsx(file)
  catalog = catalog[catalog[,sigset]==1,]

  GENELISTS <<- readRDS(sigdbgenelist)

  fcmat2 = FCMAT2_lowconc
  if(do.cov) {
    covmat = cov(fcmat2)
    COVMAT <<- covmat
  }
  covmat = COVMAT

  siglist = sort(unique(catalog$parent))
  allgenes = colnames(fcmat2)

  cat("start calculating cutoffs\n")
  if(mc.cores > 1){
    cl = makePSOCKcluster(mc.cores)
    clusterExport(cl,c("GENELISTS","COVMAT","FCMAT2_lowconc","rowMean"))
    res = parLapply(cl = cl, X=siglist, fun=cutoffCalc.inner.fc, catalog=catalog,allgenes=allgenes,pval=pval)
    cutoffs = do.call(rbind.data.frame,res)
  }
  else {
    res = lapply(X=siglist,FUN=cutoffCalc.inner.fc,catalog,allgenes,pval)
    cutoffs = do.call(rbind.data.frame,res)
  }
  if(mc.cores > 1) stopCluster(cl)

  cat("finish calculating cutoffs\n")

  if(do.compare) {
    compareSignatureCutoffs(exportSignatureCutoffs(signature_cr) ,cutoffs)
  }

  return(cutoffs)
}

#####################################################################################################
#'
#' Generates a plot comparing two cutoff calculation methods
#'
#' @param empirical_cutoffs output from the exportSignatureCutoffs function
#' @param cutoffs calculated cutoffs provided by cutoffCalc function
#' @param nsig number of signatures to compare
#' @importFrom ggplot2 ggplot
#'
#' @return ggplot object which is the scatterplot comparing the two methods
#' @export compareSignatureCutoffs

compareSignatureCutoffs = function(empirical_cutoffs, cutoffs, nsig=NULL){

  #par(mfrow=c(2,1),mar=c(5,6,6,3))
  oldcuts <- data.table(empirical_cutoffs)
  cutoffs <- data.table(cutoffs)

  #grab only signatues shared between methods
  sigs <- intersect(cutoffs$signature, oldcuts$signature)

  #sort
  cutoffs <- cutoffs[signature %in% sigs,]
  oldcuts <- oldcuts[signature %in% sigs,]

  #order
  cutoffs <- cutoffs[order(signature),]
  oldcuts <- oldcuts[order(signature),]

  dat <- data.table(empirical_method = oldcuts$cutoff, fc_method = cutoffs$cutoff)

  p <- ggplot(data = dat, aes(x = empirical_method, y = fc_method)) +
    geom_point() +
    geom_abline(slope = 1, linetype = "dashed") +
    stat_smooth(method = "lm")

  return(p)
}


#####################################################################################################
#' Inner function for the cutoff calculation based on the analytical method
#'
#' @param parent The name of the signature parent for which the cutoff is to be calculated
#' @param catalog The signature catalog
#' @param allgenes The list of all the genes in the data set
#' @param pval The p-value for the baseline distribution
#' @importFrom stats qnorm
#'
#' @return vector containing the parent (signature), cutoff, sd, bmed
#' @export cutoffCalc.inner.fc
#####################################################################################################

cutoffCalc.inner.fc = function(parent, catalog, allgenes, pval) {
  name.list = c("signature","sd","bmed","cutoff","pvalue","numsd")
  cutoffs = as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(cutoffs) = name.list
  cutoffs[1,"signature"] = parent
  cutoffs[1,"sd"] = 0
  cutoffs[1,"bmed"] = 0
  cutoffs[1,"cutoff"] = 0
  cutoffs[1,"pvalue"] = pval
  cutoffs[1,"numsd"] = 1

  scale = 1-pval/2
  IG = vector(length=length(allgenes),mode="integer")
  IG[] = 0

  parent = parent
  temp = catalog[is.element(catalog$parent,parent),]
  if(nrow(temp)==1) {
    sig = temp[1,"signature"]
    genes = GENELISTS[[sig]]
    genes.in = genes[is.element(genes,allgenes)]
    genes.out = allgenes[!is.element(allgenes,genes.in)]
    min = length(genes.in)
    mout = length(genes.out)
    IGin = IG
    IGin[is.element(allgenes,genes.in)] = 1
    IGout = IG
    IGout[is.element(allgenes,genes.out)] = 1
    wj = IGin / min - IGout / mout
    sd = sqrt(t(wj)%*%COVMAT%*%wj)
    cutoffs[1,"sd"] = sd
    cutoff = qnorm(scale,mean=0,sd=sd)
    cutoffs[1,"cutoff"] = cutoff

    mat.in = FCMAT2_lowconc[,genes.in]
    mat.out = FCMAT2_lowconc[,genes.out]
    rs.in = rowMean(mat.in)
    rs.out = rowMean(mat.out)
    diff = rs.in-rs.out
    cutoffs[1,"bmed"] = mean(diff)
  }
  else if (nrow(temp)==2) {
    sigup = temp[is.element(temp$direction,"up"),"signature"]
    upgenes = GENELISTS[[sigup]]
    upgenes.in = upgenes[is.element(upgenes,allgenes)]
    upgenes.out = allgenes[!is.element(allgenes,upgenes.in)]
    upmin = length(upgenes.in)
    upmout = length(upgenes.out)

    sigdn = temp[is.element(temp$direction,"dn"),"signature"]
    dngenes = GENELISTS[[sigdn]]
    dngenes.in = dngenes[is.element(dngenes,allgenes)]
    dngenes.out = allgenes[!is.element(allgenes,dngenes.in)]
    dnmin = length(dngenes.in)
    dnmout = length(dngenes.out)

    upIGin = IG
    upIGin[is.element(allgenes,upgenes.in)] = 1
    upIGout = IG
    upIGout[is.element(allgenes,upgenes.out)] = 1

    dnIGin = IG
    dnIGin[is.element(allgenes,dngenes.in)] = 1
    dnIGout = IG
    dnIGout[is.element(allgenes,dngenes.out)] = 1

    wj = upIGin / upmin - upIGout / upmout - dnIGin / dnmin + dnIGout / dnmout
    sd = sqrt(t(wj)%*%COVMAT%*%wj)
    cutoffs[1,"sd"] = sd

    cutoff = qnorm(scale,mean=0,sd=sd)
    cutoffs[1,"cutoff"] = cutoff

    mat.in = FCMAT2_lowconc[,upgenes.in]
    mat.out = FCMAT2_lowconc[,upgenes.out]
    uprs.in = rowMean(mat.in)
    uprs.out = rowMean(mat.out)
    mat.in = FCMAT2_lowconc[,dngenes.in]
    mat.out = FCMAT2_lowconc[,dngenes.out]
    dnrs.in = rowMean(mat.in)
    dnrs.out = rowMean(mat.out)
    diff = uprs.in-uprs.out - dnrs.in+dnrs.out
    cutoffs[1,"bmed"] = mean(diff)
  }
  return(as.vector(cutoffs))
}

#####################################################################################################
#' Utility function used by cutoffCalc.inner.fc
#'
#' @param x = matrix or data.frame
#'
#' @return vector of mean value of each row
#' @export rowMean
#####################################################################################################

rowMean <- function(x) {
  ret <- apply(x,FUN=mean,MARGIN=1)
}
