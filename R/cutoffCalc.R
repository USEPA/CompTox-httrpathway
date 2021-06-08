#' Calculate the signature-wise cutoffs based on the Setzer method
#' which does not break any correlations between genes
#'
#' @param basedir Directory that holds FCMAT2 and CHEM_DICT files.
#' @param dataset Name of actual dataset to base null data on.
#' @param sigcatalog The name of the signature catalog to use
#' @param sigset THe signature set
#' @param method The scoring method, either fc or gsea
#' @param p The p-value for the null distribution
#' @param seed Random seed.
#' @param nlowconc If not NULL, only include the lowest nlowconc concentrations for each chemical
#' @param dtxsid.exclude dtxsids to exclude, default NULL
#' @param do.load If TRUE, reload the FCMAT2 matrix, signature catalog and chemical dictionary, and store in globals
#' @param do.cov If TRUE, calculate the covariance matrix and store in a global
#' @param do.compare If TRUE, compare the cutoffs with those from the original method with no gene-gene correlation
#' @param to.file If TRUE, and do.compare=TRUE, send a plot of the comparison to a file
#' @param verbose If TRUE, write a line for each signature to show progress.
#'
#' @return No output.
#' @export
cutoffCalc = function(basedir="../input/fcdata/",
                      dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                      sigcatalog="signatureDB_master_catalog 2021-05-10",
                      sigset="screen_large",
                      method="fc",
                      pval=0.05,
                      seed = 12345,
                      nlowconc=2,
                      dtxsid.exclude=NULL,
                      do.load=T,
                      do.cov=T,
                      do.compare=F,
                      to.file=F,
                      verbose=F){
  printCurrentFunction()
  set.seed(seed)

  if(do.load) {
    #load chem_dict
    file = paste0(basedir,"CHEM_DICT_",dataset,".RData")
    load(file)
    chem_dict_filter = CHEM_DICT[CHEM_DICT$conc_index<=nlowconc,]
    if(!is.null(dtxsid.exclude)) chem_dict_filter = chem_dict_filter[!is.element(chem_dict_filter$dtxsid,dtxsid.exclude),]

    #load fcmat
    file = paste0(basedir,"FCMAT2_",dataset,".RData")
    cat("  file:",file,"\n")
    load(file)
    sk.list = chem_dict_filter$sample_key
    FCMAT2[is.nan(FCMAT2)] = 0
    FCMAT2[is.na(FCMAT2)] = 0
    FCMAT2_lowconc <<- FCMAT2[sk.list,]

    file = paste0("../input/signatures/",sigcatalog,".xlsx")
    catalog = read.xlsx(file)
    catalog = catalog[catalog[,sigset]==1,]
    file = paste0("../input/signatures/signatureDB_genelists.RData")
    print(file)
    load(file=file)
    CATALOG <<- catalog
    GENELISTS <<- genelists
  }
  catalog = CATALOG
  genelists = GENELISTS
  fcmat2 = FCMAT2_lowconc
  if(do.cov) {
    covmat = cov(fcmat2)
    COVMAT <<- covmat
  }
  covmat = COVMAT

  siglist = sort(unique(catalog$parent))
  name.list = c("signature","cutoff")
  cutoffs = as.data.frame(matrix(nrow=length(siglist),ncol=length(name.list)))
  names(cutoffs) = name.list
  cutoffs[,1] = siglist
  cutoffs[,2] = 0
  rownames(cutoffs) = siglist
  allgenes = colnames(fcmat2)

  nsig = length(siglist)
  start = 1
  #start = 16
  #nsig = 100
  scale = 0.975
  scale = 1-pval/2
  IG = vector(length=length(allgenes),mode="integer")
  IG[] = 0
  for(j in start:nsig) {
    parent = siglist[j]
    if(verbose) cat(j,nsig,parent,"\n")
    else {
      if(j%%1000==0) cat(j,nsig,parent,"\n")
    }
    temp = catalog[is.element(catalog$parent,parent),]
    if(nrow(temp)==1) {
      sig = temp[1,"signature"]
      genes = genelists[[sig]]
      genes.in = genes[is.element(genes,allgenes)]
      genes.out = allgenes[!is.element(allgenes,genes.in)]
      min = length(genes.in)
      mout = length(genes.out)
      IGin = IG
      IGin[is.element(allgenes,genes.in)] = 1
      IGout = IG
      IGout[is.element(allgenes,genes.out)] = 1
      wj = IGin / min - IGout / mout
      sd = sqrt(t(wj)%*%covmat%*%wj)
      cutoff = qnorm(scale,mean=0,sd=sd)
      cutoffs[j,"cutoff"] = cutoff
    }
    else if (nrow(temp)==2) {
      sigup = temp[is.element(temp$direction,"up"),"signature"]
      upgenes = genelists[[sigup]]
      upgenes.in = upgenes[is.element(upgenes,allgenes)]
      upgenes.out = allgenes[!is.element(allgenes,upgenes.in)]
      upmin = length(upgenes.in)
      upmout = length(upgenes.out)

      sigdn = temp[is.element(temp$direction,"dn"),"signature"]
      dngenes = genelists[[sigdn]]
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
      sd = sqrt(t(wj)%*%covmat%*%wj)

      cutoff = qnorm(scale,mean=0,sd=sd)
      cutoffs[j,"cutoff"] = cutoff
    }
    if(cutoffs[j,"cutoff"]==0) browser()
  }
  file = paste0("../output/signature_cutoff/signature_cutoff_",sigset,"_",dataset,"_",method,"_",pval,"_",nlowconc,"_with_gene_correlations.xlsx")
  write.xlsx(cutoffs,file)

  if(do.compare) {
    if(to.file) {
      file = paste0("../output/signature_cutoff/signature_cutoff_",sigset,"_",dataset,"_",method,"_",pval,"_",nlowconc,"_with_gene_correlations.pdf")
      pdf(file=file,width=5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(2,1),mar=c(5,6,6,3))

    file = paste0("../output/signature_cutoff/signature_cutoff_",sigset,"_",dataset,"_",method,"_",pval,".xlsx")
    print(file)
    oldcuts = read.xlsx(file)
    rownames(oldcuts) = oldcuts$signature
    oldcuts = oldcuts[siglist[1:nsig],]
    x = oldcuts$cutoff
    y = cutoffs[1:nsig,"cutoff"]
    plot(y~x,xlim=c(0,0.25),ylim=c(0,0.25),xlab="old method",ylab="new method",pch=".")
    lines(c(0,0.25),c(0,0.25))

    if(!to.file) browser()
    else dev.off()
  }
}
