#--------------------------------------------------------------------------------------
#'
#' Convert the signature conc resp RData file to Excel
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#'
#' u2os_toxcast_pfas_pe1_normal_refchems
#' heparg2d_toxcast_pfas_pe1_normal_refchems
#'
#' u2os_toxcast_pfas_pe1_normal_refchems
#' heparg2d_toxcast_pfas_pe1_normal_refchems
#--------------------------------------------------------------------------------------
signatureConcRespToExcel <- function(do.load=F,
                                     activeOnly=T,
                                     dataset="u2os_toxcast_pfas_pe1_normal",
                                     sigset="stress_test",
                                     sigcatalog="signatureDB_master_catalog 2021-02-17",
                                     method="fc") {
  printCurrentFunction()

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    mat = mat[mat$hitcall>0.9,]
    mat = mat[mat$top_over_cutoff>1.5,]
    SIGNATURE_CR <<- mat
  }
  mat = SIGNATURE_CR
  if(activeOnly) {
    mat = mat[mat$hitcall>0.9,]
    mat = mat[mat$top_over_cutoff>1.5,]
  }
  file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.xlsx")
  if(activeOnly)ile = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits_active_only.xlsx")
  write.xlsx(mat,file)
}
