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
#--------------------------------------------------------------------------------------
signatureChemStats <- function(do.read=F,
                               dataset="u2os_toxcast_pfas_pe1_normal",
                               sigset="screen_large",
                               method="fc",
                               celltype="U2OS") {
  printCurrentFunction(paste0(dataset,sigset,method))
  if(do.read) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR
    MAT <<- mat
  }
  mat <- MAT
  mat1 <- mat
  mat1 <- mat1[mat1$bmd<1,]
  mat1 <- mat1[mat1$hitcall>0.95,]
  mat1 <- mat1[mat1$top_over_cutoff>2,]
  sigs0 <- mat$signature
  sigs1 <- mat1$signature
  sig.list <- sort(unique(mat$signature))
  name.list <- c("celltype","signature","nall","nhigh","fhigh")
  nsig <- length(sig.list)
  res <- as.data.frame(matrix(nrow=nsig,ncol=length(name.list)))
  names(res) <- name.list
  res$celltype <- celltype
  n0 <- length(sigs0[is.element(sigs0,sig.list[1])])
  for(i in 1:nsig) {
    sig <- sig.list[i]
    n1 <- length(sigs0[is.element(sigs1,sig)])
    res[i,"signature"] <- sig
    res[i,"nall"] <- n0
    res[i,"nhigh"] <- n1
    res[i,"fhigh"] <- n1/n0
    if(i%%100==0) {
      cat("finished:",i," out of ",nsig,"\n")
    }
  }
  file <- paste0("../output/signature_cluster/",celltype,"/signature_stats_",celltype,"_",sigset,"_",dataset,"_",method,".xlsx")
  print(file)
  write.xlsx(res,file)
}
