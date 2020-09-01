#--------------------------------------------------------------------------------------
#' Generate the PFAS input data set
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' u2os_toxcast_pfas_pe1_normal
#--------------------------------------------------------------------------------------
pfasDatasetPrep <- function(do.load=F,
                             dataset="heparg2d_toxcast_pfas_pe1_normal",
                             sigset="screen_large",
                             method="fc",
                             celltype="HepaRG") {
  printCurrentFunction(paste(dataset,sigset,method))
  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }
  mat = MAT

  auc = abs(mat$top) * (3-log10(mat$bmd))
  auc[auc<0] = 0
  mat$auc = auc
  mat$ci <- exp(mat$er)*2.7765

  file <- "../input/pfas/QC sample map 2020-05-04.xlsx"
  qc <- read.xlsx(file)
  spid.list <- qc$spid
  mat <- mat[is.element(mat$sample_id,spid.list),]
  dataset = paste0("PFAS_",celltype)
  file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  SIGNATURE_CR = mat
  save(SIGNATURE_CR,file=file)

  file = paste0("../output/signature_refchemdb/validated_signatures_merged.xlsx")
  temp=read.xlsx(file)
  sig.list = temp$signature
  mat = mat[is.element(mat$signature,sig.list),]
  mat = mat[mat$hitcall>=0.95,]
  dataset = paste0("PFAS_",celltype,"_filtered")
  file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  SIGNATURE_CR = mat
  save(SIGNATURE_CR,file=file)
}

