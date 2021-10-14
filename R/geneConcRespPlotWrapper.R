#' Wrapper for all of the conc-response plotting o genes
#'
#' @param dataset Name of the data set.
#' @param mc.cores Number of cores to parallelize with.
#' @param do.load If TRUE, load the SIGNATURE_CR file, otherwiseassume that it is in memory
#' @param to.file to.file = T saves the output to a file; otherwise it's returned.
#' @param pval Desired cutoff p-value.
#' @param plotrange The x-range of the plot as a vector of 2 elements, this can be changed for special cases, but defaults to 0.001 to 100
#' @param onefile If TRUE, put all plots into one file, instead of one filer per chemical
#' @param chemfile A file of chemicals to use. If NULL, plot all chemicals
#'
#' @export
#----------------------------------------------------------------------------------
geneConcRespPlotWrapper <- function(dataset="tox21_cpp5_heparg_pe1_normal",
                                    mc.cores=20,
                                    do.load=T,
                                    to.file=F,
                                    pval = .05,
                                    plotrange=c(0.0001,100),
                                    onefile=T,
                                    chemfile=NULL) {

  printCurrentFunction(paste(dataset))
  if(do.load) {
    file <- paste0("../output/gene_conc_resp_summary/GENE_CR_",dataset,"_", pval,"_conthits.RData")
    print(file)
    load(file=file)

    GENE_CR <<- GENE_CR
  }

  #GENE_CR = GENE_CR[is.element(GENE_CR$gene,genes),]
  if(!is.null(chemfile)) {
    temp = read.xlsx(chemfile)
    dtxsid.list =temp$dtxsid
    GENE_CR = GENE_CR [is.element(GENE_CR$dtxsid,dtxsid.list),]
  }
  #fix chemical name so it can be part of a file name
  GENE_CR$proper_name = gsub("\\)","",GENE_CR$name)
  GENE_CR$proper_name = gsub("\\(","",GENE_CR$proper_name)
  GENE_CR$proper_name = gsub(":","",GENE_CR$proper_name)
  GENE_CR$proper_name = gsub("%","Percent",GENE_CR$proper_name)

  temp = GENE_CR[GENE_CR$hitcall>0.9,]
  temp = temp[temp$top_over_cutoff>1.5,]
  if(nrow(temp)<10) temp <- GENE_CR[1:10,]
  GENE_CR <- temp

  dir.create("../output/gene_conc_resp_plots/", showWarnings = F)
  foldname = paste0("../output/sgene_conc_resp_plots/",dataset,"_", pval,"_conthits")
  dir.create(foldname, showWarnings = F)
  pnames = unique(GENE_CR$proper_name)

  if(onefile) {
    if(to.file) {
      fname <- paste0("../output/gene_conc_resp_plots/gene_conc_resp_",dataset,"_",pval,"_conthits.pdf")
      pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(3,2),mar=c(5,5,4,2))
    GENE_CR = GENE_CR[order(GENE_CR$name),]
    GENE_CR = GENE_CR[order(GENE_CR$gene),]
    GENE_CR$auc <- abs(GENE_CR$top) * (3-log10(GENE_CR$bmd))
    for(i in 1:nrow(GENE_CR)) {
      geneConcRespPlot(GENE_CR[i,],plotrange=plotrange)
      if(!to.file) browser()
    }
    if(to.file) dev.off()
  }
  else {
    #cycle through chemicals for plotting (each gets its own file)
    cat("start the plotting \n")
    if(mc.cores > 1){
      cl = makePSOCKcluster(mc.cores)
      clusterExport(cl, c("plotouterGene", "geneConcRespPlot"))
      clusterExport(cl, c("acy", "acgnlsobj", "bmdbounds", "bmdobj", "cnst", "exp2", "exp3", "exp4", "exp5", "fitcnst", "fithill", "fitgnls",
                          "fitcnst", "fitpoly1", "fitpoly2", "fitpow", "fitexp2", "fitexp3","fitexp4", "fitexp5",
                          "gnls" , "gnlsderivobj", "hillfn", "hitcontinner","hitloginner","tcplfit2_core","tcplhit2_core", "loggnls", "loghill", "nestselect",
                          "poly1", "poly2", "pow", "tcplObj", "toplikelihood"))

      cat("run with cores:",mc.cores,"\n")
      output = clusterEvalQ(cl, library(stringr))
      output = parLapply(cl = cl, X=as.list(pnames), fun=plotouterGene,
                         GENE_CR = GENE_CR, foldname = foldname, plotrange=plotrange)
    } else {
      cat("run with cores:",mc.cores,"\n")
      output = lapply(X=as.list(pnames), plotouterGene,GENE_CR = GENE_CR, foldname = foldname, plotrange=plotrange)
    }

    if(mc.cores > 1) stopCluster(cl)
  }

}

