#--------------------------------------------------------------------------------------
#' Pull out the most extreme hits
#'
#' @param dataset The L2fc matrix data set
#' @param dir The directory where the data file lives
#' @param do.read If TRUE, read in FCMAT2 to a gloabal
#' @return
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_toxcast_pfas_pe1_normal
#' DMEM_6hr_pilot_normal_pe_1
#--------------------------------------------------------------------------------------
signatureChemBigHits <- function(to.file=F,
                                 do.read=F,
                                 dataset="heparg2d_toxcast_pfas_pe1_normal",
                                 sigset="screen_large",
                                 method="fc",
                                 celltype="HepaRG") {
  printCurrentFunction(paste0(dataset,sigset,method))
  if(do.read) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR
    mat <- mat[mat$hitcall>=0.95,]
    mat <- mat[mat$bmd<1,]
    mat <- mat[mat$top_over_cutoff>2,]
    MAT <<- mat
  }
  mat <- MAT
  mat$auc <- abs(mat$top) * (3-log10(mat$bmd))
  cat("Rows:",nrow(mat),"\n")

  ##################################################################################
  ##################################################################################
  ##################################################################################
  file <- paste0("../output/signature_cluster/",celltype,"/signature_bighits_",celltype,"_",sigset,"_",dataset,"_",method,".xlsx")
  write.xlsx(mat,file)

  ##################################################################################
  ##################################################################################
  ##################################################################################
  if(to.file) {
    fname <- paste0("../output/signature_cluster/",celltype,"/signature_bighits_",celltype,"_",sigset,"_",dataset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  mat$proper_name <- mat$name
  name.list <- sort(unique(mat$name))
  for(name in name.list) {
    cat(name,"\n")
    temp <- mat[is.element(mat$name,name),]
    temp <- temp[order(temp$auc,decreasing=T),]
    if(nrow(temp)>10) temp <- temp[1:10,]
    for(i in 1:nrow(temp)) {
      signatureConcRespPlot(temp[i,])
      if(!to.file) browser()
    }
  }
  if(to.file) dev.off()
}
