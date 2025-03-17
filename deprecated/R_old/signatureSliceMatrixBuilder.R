#--------------------------------------------------------------------------------------
#' Turn the signatuer slices into rectangular matrices
#'
#'   heparg2d_toxcast_pfas_pe1_normal
#'   mcf7_ph1_pe1_normal_block_123_allPG
#'   u2os_toxcast_pfas_pe1_normal
#'   PFAS_HepaRG
#'   PFAS_U2OS
#'   u2os_pilot_pe1_normal_null_pilot_lowconc
#'   u2os_toxcast_pfas_pe1_normal_refchems
#'   heparg2d_toxcast_pfas_pe1_normal_refchems
#'
#' run signatureSliceBuilder before this
#--------------------------------------------------------------------------------------
signatureSliceMatrixBuilder <- function(dataset="u2os_toxcast_pfas_pe1_normal",
                                        sigset="screen_large",
                                        method="fc",
                                        celltype="U2OS",
                                        hccut=0.95,
                                        tccut=1.5) {
  printCurrentFunction(paste(dataset,sigset,method))

  file = paste0("../output/signature_conc_resp_summary/SLICE_",sigset,"_",dataset,"_",method,"_",hccut,"_",tccut,".RData")
  load(file=file)
  mat = SLICE
  cat(nrow(mat),"\n")
 # mat = mat[1:100000,]
  mat$row_index = paste(mat$sample_id,mat$dtxsid,mat$conc)
  chems = mat[,c("row_index","sample_id","dtxsid","name","conc","celltype")]

  temp = reshape2::dcast(mat,row_index~signature,value.var="resp",fill=0,fun.aggregate=mean)
  rownames(temp) = temp[,1]
  temp = temp[,2:ncol(temp)]
  SLICEMAT = temp
  print(dim(SLICEMAT))
  SLICEMATCHEMS = chems
  file = paste0("../output/signature_conc_resp_summary/SLICEMAT_",sigset,"_",dataset,"_",method,"_",hccut,"_",tccut,".RData")
  save(SLICEMAT,file=file)
  file = paste0("../output/signature_conc_resp_summary/SLICEMATCHEMS_",sigset,"_",dataset,"_",method,"_",hccut,"_",tccut,".RData")
  save(SLICEMATCHEMS,file=file)
}

