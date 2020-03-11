library(openxlsx)
library(tcplfit2)
#' Pathway Concentration Response (P-value)
#'
#' Performs signature concentration response using p-value based cutoffs.
#'
#' Null dataset and dataset should have already been scored using signatureScore
#' and the given sigset and method. This function prepares signatureScore output
#' for CR processing, calls signatureConcRespCore_pval, formats the output,
#' saves it to disk, then calls plotouter for CR plots, if desired. If
#' conthits = T and nametag = NULL, the nametag "conthits" is automatically added
#' to the output file.
#'
#' @param sigset Name of the signature set.
#' @param dataset Name of the data set.
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#' @param nullset Name of the null data set.
#' @param mc.cores Number of cores to parallelize with.
#' @param to.file to.file = T saves the output to a file; otherwise it's returned.
#' @param do.plot do.plot = T creates concentration-response plots for every
#'   sample/signature combination and saves to disk.
#' @param pval Desired cutoff p-value.
#' @param nametag Optional descriptor tag to attach to file outputs for
#'   experimental/non-default runs.
#' @param conthits conthits = T uses continuous hitcalls, otherwise it's discrete.
#' @param aicc aicc = T uses corrected AIC to choose winning method; otherwise
#'   regular AIC.
#' @param minsigsize Minimum allowed signature size. Sample/signature combinations
#'   with less than this number of non-missing l2fc's will be discarded.
#' @param fitmodels Vector of model names to use. Probably should include "cnst".
#'
#' @import data.table
#' @import parallel
#' @import openxlsx
#'
#' @return If to.file = T, nothing. If to.file = F, dataframe with signature CR
#'   output.
#' @export
signatureConcResp <- function(sigset="pilot_tiny",
                              sigcatalog="signatureDB_master_catalog 2020-01-31",
                              dataset="DMEM_6hr_screen_normal_pe_1",
                              method="mygsea",
                              nullset="DMEM_6hr_screen_normal_pe_1_RAND1000",
                              mc.cores=1,
                              to.file=T,
                              do.plot = F,
                              pval = .05,
                              nametag = NULL,
                              conthits = T,
                              aicc = F,
                              minsigsize = 10,
                              fitmodels = c("cnst", "hill", "gnls", "poly1", "poly2", "pow", "exp2", "exp3",
                                            "exp4", "exp5"),
                              CYTOTOX) {

  printCurrentFunction(paste(dataset,sigset,method))
  starttime = proc.time()

  annotations <- signatureCatalogLoader(sigset,sigcatalog)

  #automatically add conthits nametag, then add _ to nametag if it exists
  if(is.null(nametag) && conthits) nametag = "conthits"
  if(!is.null(nametag)) nametag = paste0("_", nametag)

  # get signaturescoremat
  file <- paste0("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,".RData")
  load(file)
  #remove scores under path limit
  signaturescoremat = signaturescoremat[signaturescoremat$size >= minsigsize,]

  #attach pval cutoff and bmed
  pvalkey = getpvalcutoff(sigset, nullset, method, pvals = pval, numsds = 1)
  cat("> signatureConcResp 1\n")
  signaturescoremat$cutoff = pvalkey$cutoff[!is.na(pvalkey$pvalue)][match(signaturescoremat$signature, pvalkey$signature)]
  signaturescoremat$bmed = pvalkey$bmed[!is.na(pvalkey$pvalue)][match(signaturescoremat$signature, pvalkey$signature)]
  signaturescoremat$onesd = pvalkey$cutoff[!is.na(pvalkey$numsd)][match(signaturescoremat$signature, pvalkey$signature)]
  cat("> signatureConcResp 2\n")

  #aggregate signaturescoremat by unique sample/signature per row; data table is considerably faster than aggregate
  if(is.element("time",names(signaturescoremat)))
    signaturescoremat = setDT(signaturescoremat)[, list(conc = list(conc),resp = list(signature_score), size = min(size)),
                                     by = list(sample_id, dtxsid, casrn, name, time, signature, bmed, cutoff, onesd)]
  else
    signaturescoremat = setDT(signaturescoremat)[, list(conc = list(conc),resp = list(signature_score), size = min(size)),
                                                 by = list(sample_id, dtxsid, casrn, name, signature, bmed, cutoff, onesd)]

  cat("> signatureConcResp 3\n")

  colnames(signaturescoremat)[colnames(signaturescoremat) == "size"] = "signature_size"
  ordering = order(tolower(signaturescoremat$name), tolower(signaturescoremat$signature))
  signaturescoremat = signaturescoremat[ordering,]
  cat("> signatureConcResp 4\n")

  #turn signaturescoremat into a list of rows for lapply to use
  signaturescoremat = as.list(as.data.frame(t(signaturescoremat), stringsAsFactors = F))
  cat("> signatureConcResp 5\n")

  #lapply with inner function: signatureConcRespCore_pval
  if(mc.cores > 1){
    cat("> signatureConcResp 6\n")
    cl = makePSOCKcluster(mc.cores)
    # hideout = clusterCall(cl, function(){source("R/signatureConcRespCore_pval.R")})
    clusterExport(cl, c("acy", "acgnlsobj", "bmdbounds", "bmdobj", "cnst", "exp2", "exp3", "exp4", "exp5", "fitcnst", "fithill", "fitgnls",
                        "fitcnst", "fitpoly1", "fitpoly2", "fitpow", "fitexp2", "fitexp3","fitexp4", "fitexp5",
                        "gnls" , "gnlsderivobj", "hillfn", "hitcontinner","hitloginner","tcplfit2_core","tcplhit2_core", "loggnls", "loghill", "nestselect",
                        "poly1", "poly2", "pow", "tcplObj", "toplikelihood"))
    #SIGNATURE_CR = parLapplyLB(cl = cl, X=signaturescoremat, fun=signatureConcRespCore_pval, fitmodels = fitmodels,
    #                         conthits =conthits, aicc = aicc, chunk.size = ceiling(length(signaturescoremat)/5/mc.cores) )
    SIGNATURE_CR = parLapplyLB(cl = cl, X=signaturescoremat, fun=concRespCore, fitmodels = fitmodels,
                             conthits =conthits, aicc = aicc, verbose=FALSE, chunk.size = ceiling(length(signaturescoremat)/5/mc.cores) )
    cat("> signatureConcResp 7\n")
  } else {
    cat("> signatureConcResp 6\n")
    #SIGNATURE_CR = lapply(X=signaturescoremat, FUN = signatureConcRespCore_pval, fitmodels = fitmodels, conthits= conthits, aicc = aicc)
    SIGNATURE_CR = lapply(X=signaturescoremat, FUN =concRespCore, fitmodels = fitmodels, conthits= conthits, aicc = aicc,verbose=F)
    cat("> signatureConcResp 7\n")
  }

  #construct SIGNATURE_CR
  SIGNATURE_CR = as.data.frame(rbindlist(SIGNATURE_CR))
  rm(signaturescoremat)
  cat("> signatureConcResp 8\n")

  #add details to SIGNATURE_CR
  SIGNATURE_CR$target_class <- annotations$target_class[match(SIGNATURE_CR$signature, annotations$parent)]
  SIGNATURE_CR$super_target <- annotations$super_target[match(SIGNATURE_CR$signature, annotations$parent)]
  SIGNATURE_CR$sigset = sigset
  SIGNATURE_CR$dataset = dataset
  SIGNATURE_CR$method = method
  SIGNATURE_CR$ac50_loss = as.numeric(SIGNATURE_CR$ac50_loss)
  cat("> signatureConcResp 9\n")

  #save SIGNATURE_CR
  if(to.file){
    dir.create("../output/signature_conc_resp_summary/", showWarnings = F)
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_", pval, nametag ,".RData")
    save(SIGNATURE_CR,file=file)
  }
  cat("> signatureConcResp 10\n")

  #plotting
  if(do.plot){
    #fix chemical name so it can be part of a file name
    SIGNATURE_CR$proper_name = gsub("\\)","",SIGNATURE_CR$name)
    SIGNATURE_CR$proper_name = gsub("\\(","",SIGNATURE_CR$proper_name)
    SIGNATURE_CR$proper_name = gsub(":","",SIGNATURE_CR$proper_name)
    SIGNATURE_CR$proper_name = gsub("%","Percent",SIGNATURE_CR$proper_name)

    dir.create("../output/signature_conc_resp_plots/", showWarnings = F)
    foldname = paste0("../output/signature_conc_resp_plots/",sigset,"_",dataset,"_",method,"_", pval, nametag)
    dir.create(foldname, showWarnings = F)
    pnames = unique(SIGNATURE_CR$proper_name)

    #cycle through chemicals for plotting (each gets its own file)
    if(mc.cores > 1){
      clusterExport(cl, c("plotouter", "signatureConcRespPlot"))
      output = clusterEvalQ(cl, library(stringr))
      output = parLapply(cl = cl, X=as.list(pnames), fun=plotouter,
                         SIGNATURE_CR = SIGNATURE_CR, foldname = foldname, CYTOTOX=CYTOTOX)
    } else {
      output = lapply(X=as.list(pnames), plotouter,SIGNATURE_CR = SIGNATURE_CR, foldname = foldname, CYTOTOX=CYTOTOX)
    }

    SIGNATURE_CR$proper_name = NULL

  }
  cat("> signatureConcResp 11\n")

  if(mc.cores > 1) stopCluster(cl)
  print(proc.time() - starttime)

  if(!to.file) return(SIGNATURE_CR)
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
#' @param CYTOTOX The cytotoxicity data for all chemicals
#' @import grDevices
#'
#' @return No output.
#' @export
plotouter = function(proper_name, SIGNATURE_CR, foldname, CYTOTOX){
  #open pdf for plots
  fname <- paste0(foldname,"/conc_resp_",proper_name,".pdf")
  pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  #narrow down to given chemical
  subframe = SIGNATURE_CR[SIGNATURE_CR$proper_name == proper_name,]

  subframe = subframe[order(-subframe$hitcall, subframe$bmd),] #order by potency (optional)

  #cycle through signatures (rows) and run signatureConcRespPlot
  for(i in 1:nrow(subframe)){
    signatureConcRespPlot(subframe[i,],CYTOTOX)
  }

  graphics.off()
}

