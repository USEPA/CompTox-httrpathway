#--------------------------------------------------------------------------------------
#' find the overlap between pairs of signaturss that occur often
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
plotDEGsignatureStats <- function(to.file=F,
                                  do.read=F,
                                  dataset="u2os_toxcast_pfas_pe1_normal",
                                  sigset="screen_large",
                                  method="fc",
                                  celltype="U2OS") {
  printCurrentFunction(paste(dataset,sigset,method))
  if(do.read) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR
    #mat <- mat[mat$bmd<1,]
    mat <- mat[mat$hitcall>0.99,]
    mat <- mat[mat$top_over_cutoff>2,]
    MAT <<- mat

    file <- paste0("../output/gene_deg/",dataset,"_gene_deg.xlsx")
    DEG <<- read.xlsx(file)
  }
  if(to.file) {
    fname <- paste0("../output/gene_deg/",dataset,"_gene_signature_deg.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,2),mar=c(4,4,2,2))
  deg <- DEG
  mat <- MAT
  chem.list <- sort(unique(deg$name))
  for(chem in chem.list) {
    cat(chem,"\n")
    ctemp <- deg[is.element(deg$name,chem),]
    sid.list <- unique(ctemp$sample_id)
    for(sid in sid.list) {
      cstemp <- ctemp[is.element(ctemp$sample_id,sid),]
      res <- cstemp[,c("conc","MAD_pos4","MAD_neg4","MAD_neg5")]
      names(res) <- c("conc","DEG(pos)","DEG(neg)","BMDs")
      res$BMDs <- 0
      mtemp <- mat[is.element(mat$sample_id,sid),]
      if(nrow(mtemp)>0) {
        for(i in 1:nrow(res)) {
          c0 <- res[i,"conc"]
          if(i<nrow(res)) c1 <- res[i+1,"conc"]
          else c1 <- 1000
          m1 <- mtemp[mtemp$bmd>c0,]
          m1 <- m1[m1$bmd<=c1,]
          res[i,"BMDs"] <- nrow(m1)
        }
      }
      height <- as.matrix(res[,2:4])
      height[height==0] <- 1
      barplot(height,beside=T,cex.lab=1.2,cex.axis=1.2,main=paste0(chem,"\n",celltype," ",sid),
              ylim=c(1,10000),log="y")
      for(i in c(1,10,100,1000,10000)) lines(c(0,100),c(i,i))
      if(!to.file) browser()
    }
  }

  if(to.file) dev.off()

}
