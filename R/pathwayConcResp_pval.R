library(openxlsx)
library(tcplfit2)
#' Pathway Concentration Response (P-value)
#'
#' Performs pathway concentration response using p-value based cutoffs.
#'
#' Null dataset and dataset should have already been scored using pathwayScore
#' and the given pathset and method. This function prepares pathwayScore output
#' for CR processing, calls pathwayConcRespCore_pval, formats the output,
#' saves it to disk, then calls plotouter for CR plots, if desired. If
#' conthits = T and nametag = NULL, the nametag "conthits" is automatically added
#' to the output file.
#'
#' @param pathset Name of the pathway set.
#' @param dataset Name of the data set.
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#' @param nullset Name of the null data set.
#' @param mc.cores Number of cores to parallelize with.
#' @param to.file to.file = T saves the output to a file; otherwise it's returned.
#' @param do.plot do.plot = T creates concentration-response plots for every
#'   sample/pathway combination and saves to disk.
#' @param pval Desired cutoff p-value.
#' @param nametag Optional descriptor tag to attach to file outputs for
#'   experimental/non-default runs.
#' @param conthits conthits = T uses continuous hitcalls, otherwise it's discrete.
#' @param aicc aicc = T uses corrected AIC to choose winning method; otherwise
#'   regular AIC.
#' @param minpathsize Minimum allowed pathway size. Sample/pathway combinations
#'   with less than this number of non-missing l2fc's will be discarded.
#' @param fitmodels Vector of model names to use. Probably should include "cnst".
#'
#' @import data.table
#' @import parallel
#' @import openxlsx
#'
#' @return If to.file = T, nothing. If to.file = F, dataframe with pathway CR
#'   output.
#' @export
pathwayConcResp_pval <- function(pathset="PathwaySet_20191107",
                                 dataset="DMEM_6hr_pilot_normal_pe_1",
                                 method="mygsea",
                                 nullset = "DMEM_6hr_pilot_normal_pe_1_RAND1000",
                                 mc.cores=1,
                                 to.file=T,
                                 do.plot = F,
                                 pval = .05,
                                 nametag = NULL,
                                 conthits = T,
                                 aicc = F,
                                 minpathsize = 10,
                                 fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                               "exp4", "exp5"),
                                 pathway_annotations,
                                 CYTOTOX) {

  printCurrentFunction(paste(dataset,pathset,method))
  starttime = proc.time()

  #automatically add conthits nametag, then add _ to nametag if it exists
  if(is.null(nametag) && conthits) nametag = "conthits"
  if(!is.null(nametag)) nametag = paste0("_", nametag)

  # get pathscoremat
  file <- paste0("../output/pathway_score_summary/PATHSCOREMAT_",pathset,"_",dataset,"_",method,".RData")
  load(file)
  #remove scores under path limit
  pathscoremat = pathscoremat[pathscoremat$size >= minpathsize,]

  #attach pval cutoff and bmed
  pvalkey = getpvalcutoff(pathset, nullset, method, pvals = pval, numsds = 1)
  cat(">pathwayConcResp_pval 1\n")
  pathscoremat$cutoff = pvalkey$cutoff[!is.na(pvalkey$pvalue)][match(pathscoremat$pathway, pvalkey$pathway)]
  pathscoremat$bmed = pvalkey$bmed[!is.na(pvalkey$pvalue)][match(pathscoremat$pathway, pvalkey$pathway)]
  pathscoremat$onesd = pvalkey$cutoff[!is.na(pvalkey$numsd)][match(pathscoremat$pathway, pvalkey$pathway)]
  cat(">pathwayConcResp_pval 2\n")

  #aggregate pathscoremat by unique sample/pathway per row; data table is considerably faster than aggregate
  pathscoremat = setDT(pathscoremat)[, list(conc = list(conc),resp = list(pathway_score), size = min(size)),
                                     by = list(sample_id, dtxsid, casrn, name, time, pathway, bmed, cutoff, onesd)]
  cat(">pathwayConcResp_pval 3\n")

  colnames(pathscoremat)[colnames(pathscoremat) == "size"] = "pathway_size"
  ordering = order(tolower(pathscoremat$name), tolower(pathscoremat$pathway))
  pathscoremat = pathscoremat[ordering,]
  cat(">pathwayConcResp_pval 4\n")

  #turn pathscoremat into a list of rows for lapply to use
  pathscoremat = as.list(as.data.frame(t(pathscoremat), stringsAsFactors = F))
  cat(">pathwayConcResp_pval 5\n")

  #lapply with inner function: pathwayConcRespCore_pval
  if(mc.cores > 1){
    cl = makePSOCKcluster(mc.cores)
    # hideout = clusterCall(cl, function(){source("R/pathwayConcRespCore_pval.R")})
    clusterExport(cl, c("acy", "acgnlsobj", "bmdbounds", "bmdobj", "cnst", "exp2", "exp3", "exp4", "exp5", "fitcnst", "fithill", "fitgnls",
                        "fitcnst", "fitpoly1", "fitpoly2", "fitpow", "fitexp2", "fitexp3","fitexp4", "fitexp5",
                        "gnls" , "gnlsderivobj", "hillfn", "hitcontinner","hitloginner","tcplfit2_core","tcplhit2_core", "loggnls", "loghill", "nestselect",
                        "poly1", "poly2", "pow", "tcplObj", "toplikelihood"))
    #PATHWAY_CR = parLapplyLB(cl = cl, X=pathscoremat, fun=pathwayConcRespCore_pval, fitmodels = fitmodels,
    #                         conthits =conthits, aicc = aicc, chunk.size = ceiling(length(pathscoremat)/5/mc.cores) )
    PATHWAY_CR = parLapplyLB(cl = cl, X=pathscoremat, fun=concRespCore, fitmodels = fitmodels,
                             conthits =conthits, aicc = aicc, verbose=FALSE, chunk.size = ceiling(length(pathscoremat)/5/mc.cores) )
  } else {
    cat(">pathwayConcResp_pval 6\n")

    #PATHWAY_CR = lapply(X=pathscoremat, FUN = pathwayConcRespCore_pval, fitmodels = fitmodels, conthits= conthits, aicc = aicc)
    PATHWAY_CR = lapply(X=pathscoremat, FUN =concRespCore, fitmodels = fitmodels, conthits= conthits, aicc = aicc,verbose=T)
    cat(">pathwayConcResp_pval 7\n")

  }

  #construct PATHWAY_CR
  PATHWAY_CR = as.data.frame(rbindlist(PATHWAY_CR))
  rm(pathscoremat)
  cat(">pathwayConcResp_pval 8\n")

  #add details to PATHWAY_CR
  PATHWAY_CR$pathway_class <- pathway_annotations$super_class[match(PATHWAY_CR$pathway, pathway_annotations$pathway)]
  PATHWAY_CR$pathset = pathway_annotations$pathset[match(PATHWAY_CR$pathway, pathway_annotations$pathway)]
  PATHWAY_CR$dataset = rep(dataset, nrow(PATHWAY_CR))
  PATHWAY_CR$method = rep(method, nrow(PATHWAY_CR))
  PATHWAY_CR$ac50_loss = as.numeric(PATHWAY_CR$ac50_loss)
  cat(">pathwayConcResp_pval 9\n")

  #save PATHWAY_CR
  if(to.file){
    dir.create("../output/pathway_conc_resp_summary/", showWarnings = F)
    file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_", pval, nametag ,".RData")
    save(PATHWAY_CR,file=file)
  }
  cat(">pathwayConcResp_pval 10\n")

  #plotting
  if(do.plot){
    #fix chemical name so it can be part of a file name
    PATHWAY_CR$proper_name = gsub("\\)","",PATHWAY_CR$name)
    PATHWAY_CR$proper_name = gsub("\\(","",PATHWAY_CR$proper_name)
    PATHWAY_CR$proper_name = gsub(":","",PATHWAY_CR$proper_name)
    PATHWAY_CR$proper_name = gsub("%","Percent",PATHWAY_CR$proper_name)

    dir.create("../output/pathway_conc_resp_plots/", showWarnings = F)
    foldname = paste0("../output/pathway_conc_resp_plots/",pathset,"_",dataset,"_",method,"_", pval, nametag)
    dir.create(foldname, showWarnings = F)
    pnames = unique(PATHWAY_CR$proper_name)

    #cycle through chemicals for plotting (each gets its own file)
    if(mc.cores > 1){
      clusterExport(cl, c("plotouter", "pathwayConcRespPlot"))
      output = clusterEvalQ(cl, library(stringr))
      output = parLapply(cl = cl, X=as.list(pnames), fun=plotouter,
                         PATHWAY_CR = PATHWAY_CR, foldname = foldname, CYTOTOX=CYTOTOX)
    } else {
      output = lapply(X=as.list(pnames), plotouter,PATHWAY_CR = PATHWAY_CR, foldname = foldname, CYTOTOX=CYTOTOX)
    }

    PATHWAY_CR$proper_name = NULL

  }
  cat(">pathwayConcResp_pval 11\n")

  if(mc.cores > 1) stopCluster(cl)
  print(proc.time() - starttime)

  if(!to.file) return(PATHWAY_CR)
}

#' Plot Outer
#'
#' Calls pathwayConcResp plotting function.
#'
#' Calls pathwayConcResp plotting function for one chemical and every pathway.
#' Saves a single pdf to disk for the given chemical containing every pathway
#' CR plot.
#'
#' @param proper_name Chemical name to be used in file name.
#' @param PATHWAY_CR Dataframe output of pathwayConcResp_pval.
#' @param foldname Folder name for output file.
#' @param CYTOTOX The cytotoxicity data for all chemicals
#' @import grDevices
#'
#' @return No output.
#' @export
plotouter = function(proper_name, PATHWAY_CR, foldname, CYTOTOX){
  #open pdf for plots
  fname <- paste0(foldname,"/conc_resp_",proper_name,".pdf")
  pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  #narrow down to given chemical
  subframe = PATHWAY_CR[PATHWAY_CR$proper_name == proper_name,]

  subframe = subframe[order(-subframe$hitcall, subframe$bmd),] #order by potency (optional)

  #cycle through pathways (rows) and run pathwayConcRespPlot
  for(i in 1:nrow(subframe)){
    pathwayConcRespPlot(subframe[i,],CYTOTOX)
  }

  graphics.off()
}

