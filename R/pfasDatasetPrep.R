#--------------------------------------------------------------------------------------
#' Generate the PFAS input data set
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' u2os_toxcast_pfas_pe1_normal
#--------------------------------------------------------------------------------------
pfasDatasetPrep <- function(do.load=F,
                            dataset="heparg2d_toxcast_pfas_pe1_normal",
                            sigset="dorothea",
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

  file <- "../input/pfas/PFAS wide sample qc.xlsx"
  qc <- read.xlsx(file)
  spid.list.0 <- unique(qc$spid)
  #mat <- mat[is.element(mat$sample_id,spid.list),]

  dtxsid.list = unique(qc$dtxsid)
  mat <- mat[is.element(mat$dtxsid,dtxsid.list),]
  spid.list.1 <- unique(mat$sample_id)
  spid.missing = spid.list.1[!is.element(spid.list.1,spid.list.0)]
  cat("missing spids in QC:",spid.missing,"\n")
  temp = unique(mat[is.element(mat$sample_id,spid.missing),c("dtxsid","name","sample_id")])
  file = paste0("../input/missing_spid_in_qc ",dataset,".xlsx")
  write.xlsx(temp,file)
  dataset = paste0("PFAS_",celltype)
  file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  SIGNATURE_CR = mat
  save(SIGNATURE_CR,file=file)
}

