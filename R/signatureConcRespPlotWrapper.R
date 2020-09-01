library(openxlsx)
#' Wrapper for all of the conc-response plotting
#'

#' @param sigset Name of the signature set.
#' @param dataset Name of the data set.
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#' @param mc.cores Number of cores to parallelize with.
#' @param do.load If TRUE, load the SIGNATURE_CR file, otherwiseassume that it is in memory
#' to.file to.file = T saves the output to a file; otherwise it's returned.
#' @param pval Desired cutoff p-value.
#' @param nametag Optional descriptor tag to attach to file outputs for
#'   experimental/non-default runs.
#
#' @import data.table
#' @import parallel
#' @import openxlsx
#'
#' @export
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_toxcast_pfas_pe1_normal
#'----------------------------------------------------------------------------------
signatureConcRespPlotWrapper <- function(sigset="screen_large",
                                         dataset="u2os_toxcast_pfas_pe1_normal",
                                         sigcatalog,
                                         method="fc",
                                         mc.cores=20,
                                         do.load=T,
                                         pval = .05,
                                         nametag = NULL) {

  printCurrentFunction(paste(dataset,sigset,method))
  if(do.load) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_", pval, nametag ,"_conthits.RData")
    print(file)
    load(file=file)
    SIGNATURE_CR <<- SIGNATURE_CR
  }
  #annotations <- signatureCatalogLoader(sigset,sigcatalog)
  #SIGNATURE_CR$super_target <- annotations$super_target[match(SIGNATURE_CR$signature, annotations$parent)]
  #file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_", pval, nametag ,"_conthits.RData")
  #save(SIGNATURE_CR,file=file)
  #browser()
  #fix chemical name so it can be part of a file name
  SIGNATURE_CR$proper_name = gsub("\\)","",SIGNATURE_CR$name)
  SIGNATURE_CR$proper_name = gsub("\\(","",SIGNATURE_CR$proper_name)
  SIGNATURE_CR$proper_name = gsub(":","",SIGNATURE_CR$proper_name)
  SIGNATURE_CR$proper_name = gsub("%","Percent",SIGNATURE_CR$proper_name)

  temp <- SIGNATURE_CR[SIGNATURE_CR$hitcall>0.8,]
  if(nrow(temp)<10) temp <- SIGNATURE_CR[1:10,]
  SIGNATURE_CR <- temp

  dir.create("../output/signature_conc_resp_plots/", showWarnings = F)
  foldname = paste0("../output/signature_conc_resp_plots/",sigset,"_",dataset,"_",method,"_", pval, nametag)
  dir.create(foldname, showWarnings = F)
  pnames = unique(SIGNATURE_CR$proper_name)

  #cycle through chemicals for plotting (each gets its own file)
  cat("start the plotting \n")
  if(mc.cores > 1){
    cl = makePSOCKcluster(mc.cores)
    clusterExport(cl, c("plotouter", "signatureConcRespPlot"))
    clusterExport(cl, c("acy", "acgnlsobj", "bmdbounds", "bmdobj", "cnst", "exp2", "exp3", "exp4", "exp5", "fitcnst", "fithill", "fitgnls",
                        "fitcnst", "fitpoly1", "fitpoly2", "fitpow", "fitexp2", "fitexp3","fitexp4", "fitexp5",
                        "gnls" , "gnlsderivobj", "hillfn", "hitcontinner","hitloginner","tcplfit2_core","tcplhit2_core", "loggnls", "loghill", "nestselect",
                        "poly1", "poly2", "pow", "tcplObj", "toplikelihood"))

    cat("run with cores:",mc.cores,"\n")
    output = clusterEvalQ(cl, library(stringr))
    output = parLapply(cl = cl, X=as.list(pnames), fun=plotouter,
                       SIGNATURE_CR = SIGNATURE_CR, foldname = foldname)
  } else {
    cat("run with cores:",mc.cores,"\n")
    output = lapply(X=as.list(pnames), plotouter,SIGNATURE_CR = SIGNATURE_CR, foldname = foldname)
  }

  if(mc.cores > 1) stopCluster(cl)
}

#' Plot Outer
#'
#' Calls signatureConcResp plotting function.
#'
#' Calls signatureConcResp plotting function for one chemical and every signature.
#' Saves a single pdf to disk for the given chemical containing every signature
#' CR plot.
#'
#' @param proper_name Chemical name to be used in file name.
#' @param SIGNATURE_CR Dataframe output of signatureConcResp_pval.
#' @param foldname Folder name for output file.
#' @import grDevices
#'
#' @return No output.
#' @export
plotouter = function(proper_name, SIGNATURE_CR, foldname){
  #open pdf for plots
  #printCurrentFunction(proper_name)
  temp <- SIGNATURE_CR[SIGNATURE_CR$proper_name == proper_name,]
  sid.list <- unique(temp$sample_id)
  for(sid in sid.list) {
    fname <- paste0(foldname,"/conc_resp_",proper_name," ",sid,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    par(mfrow=c(3,2),mar=c(4,4,2,2))

    #narrow down to given chemical / and sample_id
    subframe = temp[is.element(temp$sample_id,sid),]

    subframe$auc <- abs(subframe$top) * (3-log10(subframe$bmd))
    subframe <- subframe[order(subframe$auc,decreasing=T),]

    #subframe = subframe[order(-subframe$hitcall,-subframe$top_over_cutoff, subframe$bmd),] #order by potency (optional)

    #cycle through signatures (rows) and run signatureConcRespPlot
    for(i in 1:nrow(subframe)){
      signatureConcRespPlot(subframe[i,])
    }
    dev.off()
  }
  graphics.off()
}

