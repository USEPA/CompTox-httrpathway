library(openxlsx)
#' Filter data for the reference chemicals
#'
#' @param sigset Name of the signature set.
#' @param dataset Name of the data set.
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#' @param bmr_scale	bmr scaling factor. Default = 1.349
#' @param mc.cores Number of cores to parallelize with.
#' @param do.load If TRUE, load the SIGNATURE_CR file, otherwiseassume that it is in memory
#' to.file to.file = T saves the output to a file; otherwise it's returned.
#' @param pval Desired cutoff p-value.
#' @param nametag Optional descriptor tag to attach to file outputs for
#'   experimental/non-default runs.
#'
#' @export
#' PFAS_HepaRG
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#' PFAS_HepaRG
#' PFAS_U2OS
#'----------------------------------------------------------------------------------
signatureConcRespRefchemFilter <- function(sigset="screen_large",
                                           dataset="u2os_toxcast_pfas_pe1_normal",
                                           sigcatalog,
                                           method="fc",
                                           mc.cores=20,
                                           do.load=T,
                                           pval = .05,
                                           nametag = "_conthits") {

  printCurrentFunction(paste(dataset,sigset,method))
  if(do.load) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_", pval,nametag,".RData")
    print(file)
    load(file=file)
    SIGNATURE_CR <<- SIGNATURE_CR
  }
  SIGNATURE_CR0 = SIGNATURE_CR
  file = "../input/chemicals/refchem_super_target_map.xlsx"
  refchems = read.xlsx(file)
  refchems = refchems[!is.element(refchems$stlist,"-"),]
  mat = SIGNATURE_CR[is.element(SIGNATURE_CR$dtxsid,refchems$dtxsid),]
  mat = mat[mat$hitcall>0,]
  mat = mat[mat$top_over_cutoff>1,]
  SIGNATURE_CR = mat
  file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_refchems_",method,"_", pval,nametag,".RData")
  save(SIGNATURE_CR,file=file)
  SIGNATURE_CR <<- SIGNATURE_CR0
}

