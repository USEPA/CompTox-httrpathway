#--------------------------------------------------------------------------------------
#' Output a subset of the signature conc resp data set for a subset of chemicals
#'
#'  heparg2d_toxcast_pfas_pe1_normal
#'  mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#'  PFAS_HepaRG
#'  PFAS_U2OS
#'  u2os_pilot_pe1_normal_null_pilot_lowconc
#'  u2os_toxcast_pfas_pe1_normal_refchems
#'  heparg2d_toxcast_pfas_pe1_normal_refchems
#'

#--------------------------------------------------------------------------------------
signatureConcRespSubset <- function(do.load=F,
                                    dataset="mcf7_ph1_pe1_normal_block_123",
                                    dataset.new = "MCF7_HTTr_Water_Samples",
                                    sigset="screen_large",
                                    method="fc",
                                    cnames=c(
                                      "MCF7_HTTr_Water_Sample_1",
                                      "MCF7_HTTr_Water_Sample_2",
                                      "MCF7_HTTr_Water_Sample_3",
                                      "MCF7_HTTr_Water_Sample_4",
                                      "MCF7_HTTr_Water_Sample_5",
                                      "MCF7_HTTr_Water_Sample_6",
                                      "MCF7_HTTr_Water_Sample_7",
                                      "MCF7_HTTr_Water_Sample_8",
                                      "MCF7_HTTr_Water_Sample_9",
                                      "MCF7_HTTr_Water_Sample_10",
                                      "MCF7_HTTr_Water_Sample_11",
                                      "MCF7_HTTr_Water_Sample_12",
                                      "MCF7_HTTr_Water_Sample_13",
                                      "MCF7_HTTr_Water_Sample_14",
                                      "MCF7_HTTr_Water_Sample_15",
                                      "MCF7_HTTr_Water_Sample_16",
                                      "MCF7_HTTr_Water_Sample_17",
                                      "MCF7_HTTr_Water_Sample_18"
                                    )) {
  printCurrentFunction(paste(dataset,sigset,method))

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }

  mat = MAT
  mat1 = mat[is.element(mat$name,cnames),]
  SIGNATURE_CR = mat1
  file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset.new,"_",method,"_0.05_conthits.RData")
  save(SIGNATURE_CR,file=file)
}

