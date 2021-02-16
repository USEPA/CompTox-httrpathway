#--------------------------------------------------------------------------------------
#' pull our HTTr data for the Myopathy chemicals
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
myopathy.httr <- function(do.load=F,
                          dataset="mcf7_ph1_pe1_normal_block_123",
                          sigset="screen_large",
                          method="fc",
                          celltype="MCF7",
                          hccut=0.9,
                          tccut=1.5) {
  printCurrentFunction(paste(dataset,sigset,method))

  file = "../myopathy/myopathy chemicals.xlsx"
  chems = read.xlsx(file)
  dtxsid.list = chems$dtxsid

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }

  res = MAT
  res = res[res$top_over_cutoff>tccut,]
  res = res[res$hitcall>hccut,]
  res= res[is.element(res$dtxsid,dtxsid.list),]


  name.list = c("sample_id","dtxsid","name","super_target","signature","top_over_cutoff","bmd","hitcall","top")
  res = res[,name.list]
  res$celltype = celltype
  file = paste0("../myopathy/myopathy ",celltype,".xlsx")
  write.xlsx(res,file)
}

