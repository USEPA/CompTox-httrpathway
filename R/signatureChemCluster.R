#--------------------------------------------------------------------------------------
#' Cluster the signatures and the chemicals
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
signatureChemCluster <- function(to.file=F,
                                 do.read=F,
                                 do.prepres=F,
                                 dataset="u2os_toxcast_pfas_pe1_normal",
                                 sigcatalog="signatureDB_master_catalog 2020-07-10",
                                 sigset="screen_large",
                                 method="fc",
                                 celltype="U2OS",
                                 tccut=2,
                                 hccut=0.9) {
  printCurrentFunction(paste0(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/signature_cluster/signature_cluster_",dataset,"_",sigset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  if(do.read) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    MAT <<-SIGNATURE_CR
    mat <- MAT
    mat <- mat[mat$hitcall>=0.95,]
    mat <- mat[mat$bmd<1,]
    mat <- mat[mat$top_over_cutoff>1.5,]
    starget <- sort(unique(mat$super_target))
    file <- paste0("../output/signature_cluster/extreme_super_target_",dataset,"_",sigset,"_",method,".xlsx")
    write.xlsx(starget,file)
  }
  if(!to.file) browser()
  if(do.prepres) {
    mat <- MAT
    mat <- mat[mat$hitcall>hccut,]
    mat <- mat[mat$top_over_cutoff>tccut,]
    mat <- mat[mat$bmd<10,]
    bmd <- 3-log10(mat$bmd)
    mat$auc <- mat$hitcall
    cat("dim(mat):",dim(mat),"\n")
    cat("number of unique signagtures: ",length(unique(mat$signature)),"\n")
    cat("number of unique chemicals: ",length(unique(mat$name)),"\n")
    hist(mat$auc,xlab="AUC",cex.axis=1.2,cex.lab=1.2,main=celltype)

    res <- reshape2::dcast(mat,name~signature,fill=0,value.var="auc",fun.aggregate=mean)
    rownames(res) <- res[,1]
    res <- as.matrix(res[,2:ncol(res)])
    RES <<- res
  }
  res <- RES
  temp <- res
  temp[temp!=0] <- 1
  cs <- colSums(temp)
  rs <- rowSums(temp)
  hist(cs,xlab="count by signature",cex.axis=1.2,cex.lab=1.2,main=celltype)
  hist(rs,xlab="count by chemical",cex.axis=1.2,cex.lab=1.2,main=celltype)
  ccut <- 5
  res <- res[rs>=ccut,]
  res <- res[,rs>=ccut]
  cat("dim(res):",dim(res),'\n')
  dist.sig <- dist(t(res))
  dist.chem <- dist(res)
  hist(dist.sig,xlab="signature distance matrix",cex.axis=1.2,cex.lab=1.2,main=celltype)
  hist(dist.chem,xlab="chemical distance matrix",cex.axis=1.2,cex.lab=1.2,main=celltype)

  sig.tree <- hclust(dist.sig,method="ward.D")
  chem.tree <- hclust(dist.chem,method="ward.D")

  par(mfrow=c(2,1),mar=c(4,4,2,2))
  plot(chem.tree, xlab="", sub="",main = "Chemical Clusters",labels = FALSE, hang = 0.04)
  plot(sig.tree, xlab="", sub="",main = "Signature Clusters",labels = FALSE, hang = 0.04)

  nchem <- nrow(res)
  nsig <- ncol(res)
  avgsize <- 5

  chem.cluster <- as.data.frame(cutree(chem.tree,k=round(nchem/avgsize)))
  chem.cluster <- as.data.frame(cbind(rownames(chem.cluster),chem.cluster[,1]))
  chem.cluster[,2] <- as.integer(chem.cluster[,2])
  names(chem.cluster) <- c("chemical","cluster")
  file <- paste0("../output/signature_cluster/signature_cluster_",dataset,"_",sigset,"_",method,"_chems.xlsx")
  write.xlsx(chem.cluster,file)

  sig.cluster <- as.data.frame(cutree(sig.tree,k=round(nsig/avgsize)))
  sig.cluster <- as.data.frame(cbind(rownames(sig.cluster),sig.cluster[,1]))
  sig.cluster[,2] <- as.integer(sig.cluster[,2])
  names(sig.cluster) <- c("signature","cluster")
  file <- paste0("../output/signature_cluster/signature_cluster_",dataset,"_",sigset,"_",method,"_signature.xlsx")
  write.xlsx(sig.cluster,file)


  if(!to.file) browser()
  else dev.off()
}
