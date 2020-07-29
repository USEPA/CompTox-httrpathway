library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
library(gplots)
library(RColorBrewer)
#--------------------------------------------------------------------------------------
#' Tlook at trends in the duplicated chemicals
#' @param method signature scoring method in c("fc", "gsva", "mygsea")
#'
#'
#' Error bars are exp(er)*qt(.025,4) = exp(er)*2.7765
#--------------------------------------------------------------------------------------
dup.trends <- function(to.file=F,
                       hccut=0.9,
                       tccut=2,
                       dbmdcut=30,
                       do.err.filter=T,
                       dataset="meanncnt0_5-plateteffect_1-shrinkage_normal",
                       sigcatalog="signatureDB_master_catalog 2020-05-05",
                       sigset="screen_large",
                       method="mygsea",
                       do.load=F) {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/ dup.trend ",dataset,"_",sigset," ",do.err.filter," ",hccut," ",tccut,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  if(do.load) {
    annotations <- signatureCatalogLoader(sigset,sigcatalog)
    ANN <<- annotations

    file <- paste0("../input/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits active.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR_ACTIVE
    #browser()
    mat <- mat[mat$hitcall>0.5,]
    MAT <<- mat
  }
  mat <- MAT
  mat[is.na(mat$super_target),"super_target"] <- "-"
  mat <- mat[!is.element(mat$super_target,"-"),]
  catalog <- ANN


  res1 <- unique(mat[,c("sample_id","dtxsid","casrn","name")])
  res1 <- res1[order(res1$dtxsid),]
  rownames(res1) <- res1$sample_id
  dtxsid.list <- res1$dtxsid
  dups <- duplicated(dtxsid.list)
  temp <- dtxsid.list[dups]
  res1$duplicate <- 0
  res1[is.element(res1$dtxsid,temp),"duplicate"] <- 1
  dups <- res1[res1$duplicate==1,]
  mat <- mat[is.element(mat$sample_id,dups$sample_id),]

  # filters
  hc <- mat$hitcall
  tc <- mat$top_over_cutoff
  bmd <- mat$bmd
  er <- mat$er
  delta.bmd <- mat$bmdu-mat$bmdl
  delta.bmd[!is.na(delta.bmd)] <- 100
  err.range <- exp(er)*2.7765
  err_over_cutoff <- err.range / mat$cutoff

  err.filter <- err_over_cutoff / tc
  mask1 <- err.filter
  mask1[] <- 1
  if(do.err.filter) mask1[err.filter>1] <- 0
  cat("=============================================================\n")
  cat("errcut:",format(sum(mask1)/length(mask1),digits=3),"\n")

  mask2 <- mask1
  mask2[] <- 1
  mask2[tc<tccut] <- 0
  cat("tccut:",format(sum(mask2)/length(mask2),digits=3),"\n")

  mask3 <- mask1
  mask3[] <- 1
  mask3[hc<hccut] <- 0
  cat("hccut:",format(sum(mask3)/length(mask3),digits=3),"\n")

  mask4 <- mask1
  mask4[] <- 1
  mask4[delta.bmd<dbmdcut] <- 0
  cat("bmdcut:",format(sum(mask4)/length(mask4),digits=3),"\n")

  mask <- mask1*mask2*mask3*mask4
  mat <- mat[mask==1,]
  cat("allcut:",format(sum(mask)/length(mask),digits=3),"\n")

  res <- dcast(mat,sample_id~signature,value.var="bmd",fill=1000)
  rownames(res) <- res[,1]
  res <- res[,2:ncol(res)]
  names <- dups[rownames(res),"name"]
  lres <- -log10(res/1000)

  sid.list <- rownames(lres)
  ldups <- dups[sid.list,]
  ldups <- ldups[order(ldups$name),]
  lres <- lres[ldups$sample_id,]
  names <- ldups$name
  result <- heatmap.2(as.matrix(lres),
                      margins=c(10,10),
                      dendrogram="col",
                      scale="none",
                      xlab="",
                      ylab="",
                      cexCol=0.1,
                      cexRow=0.4,
                      Rowv=F,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="target count",
                      cex.main=1,
                      labRow=names,
                      main=paste("PFAS Dups ",do.err.filter,hccut,tccut))


  if(!to.file) browser()
  else dev.off()

}
