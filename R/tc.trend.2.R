library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
library(gplots)
library(RColorBrewer)
#--------------------------------------------------------------------------------------
#' Calulate trends in conc-response data as function of top/cutoff (tc)
#' @param method signature scoring method in c("fc", "gsva", "mygsea")
#--------------------------------------------------------------------------------------
tc.trend.2 <- function(to.file=F,
                       hccut=0.9,
                       dataset="meanncnt0_5-plateteffect_1-shrinkage_normal",
                       sigcatalog="signatureDB_master_catalog 2020-05-05",
                       sigset="screen_large",
                       method="mygsea",
                       do.load=F) {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/ tc.trend.2 ",dataset,"_",sigset,".pdf")
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
  mat <- mat[mat$hitcall>hccut,]
  mat[is.na(mat$super_target),"super_target"] <- "-"
  catalog <- ANN
  file <- "../input/QC sample map 2020-05-04.xlsx"
  qc <- read.xlsx(file)
  qc <- qc[is.element(qc$score,c("H","M","No score","TBD")),]
  spid.list <- qc$spid
  mat <- mat[is.element(mat$sample_id,spid.list),]

  ##################################################################################
  file <- "../input/PFAS categories.xlsx"
  cats <- read.xlsx(file)
  #cats <- cats[!is.element(cats$category,"TxP_PFAS_perFhexyl"),]
  mat$category <- "nocat"
  dtxsid.list <- unique(mat$dtxsid)
  for(dtxsid in dtxsid.list) {
    if(is.element(dtxsid,cats$dtxsid)) {
      categ <- cats[is.element(cats$dtxsid,dtxsid),"category"]
      mat[is.element(mat$dtxsid,dtxsid),"category"] <- categ[1]
    }
  }

  hc <- mat$hitcall
  tc <- mat$top_over_cutoff
  bmd <- mat$bmd
  er <- mat$er

  hist(hc)
  hist(tc)
  hist(tc)
  plot(density(bmd),log="x")

  par(mfrow=c(3,1),mar=c(4,4,2,2))
  vals <- c(1e-3,1e-2,1e-1,1,10,100)
  nval <- length(vals)
  x <- NULL
  y <- NULL
  for(val in vals) {
    bottom <- val
    top <- val*10
    temp <- mat[mat$bmd>=bottom,]
    temp <- temp[temp$bmd<top,]
    yy <- temp$top_over_cutoff
    xx <- vector(mode="character",length=length(yy))
    xx[] <- paste0(bottom,"-",top)
    y <- c(y,yy)
    x <- c(x,xx)
  }
  boxplot(y~x,cex.lab=1.2,cex.axis=1.2,ylab="top / cutoff",xlab="BMD",ylim=c(1,6))


  vals <- c(1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)
  nval <- length(vals)
  x <- NULL
  y <- NULL
  for(val in vals) {
    bottom <- val
    top <- val+0.5
    temp <- mat[mat$top_over_cutoff>=bottom,]
    temp <- temp[temp$top_over_cutoff<top,]
    yy <- temp$bmd
    xx <- vector(mode="character",length=length(yy))
    xx[] <- paste0(bottom,"-",top)
    y <- c(y,yy)
    x <- c(x,xx)
  }
  boxplot(y~x,log="y",cex.lab=1.2,cex.axis=1.2,xlab="top / cutoff",ylab="BMD")

  name.list <- c("cutoff","signatures","targets")
  res <- data.frame(matrix(nrow=nval,ncol=length(name.list)))
  names(res) <- name.list
  res[,1] <- vals
  for(i in 1:length(vals)) {
    tccut <- vals[i]
    mat2 <- mat[mat$top_over_cutoff>tccut,]
    hc <- mat2$hitcall
    tc <- mat2$top_over_cutoff
    bmd <- mat2$bmd
    mat2 <- mat2[mat2$bmd<10,]
    n1 <- length(unique(mat2$signature))
    n2 <- length(unique(mat2$super_target))
    res[i,"signatures"] <- n1
    res[i,"targets"] <- n2
  }

  file <- paste0("../output/ tc.trend.2 ",dataset,"_",sigset,".xlsx")
  write.xlsx(res,file)

  tccut <- 2
  ##################################################################################
  mat2 <- mat[mat$top_over_cutoff>tccut,]
  mat2 <- mat2[!is.element(mat2$super_target,"-"),]
  mat2 <- mat2[!is.element(mat2$super_target,"unknown"),]
  hc <- mat2$hitcall
  tc <- mat2$top_over_cutoff
  bmd <- mat2$bmd
  #mat2 <- mat2[mat2$bmd<10,]
  ntar <- length(unique(mat2$super_target))
  nchem <- length(unique(mat2$name))
  chems <- unique(mat2$name)
  res <- as.data.frame(matrix(nrow=nchem,ncol=ntar))
  rownames(res) <- unique(mat2$name)
  names(res) <- unique(mat2$super_target)
  res[] <- 0
  for(i in 1:nchem) {
    chem <- chems[i]
    stlist <- mat2[is.element(mat2$name,chem),"super_target"]
    for(st in unique(stlist)) {
      nst <- length(stlist[is.element(stlist,st)])
      res[chem,st] <- nst
    }
  }

  #res <- res[1:10,1:10]
  res[res>5] <- 5
  rs <- rowSums(res)
  res <- res[rs>0,]
  cs <- colSums(res)
  res <- res[,cs>0]
  result <- heatmap.2(as.matrix(res),
                      margins=c(10,10),
                      dendrogram="both",
                      scale="none",
                      main="",
                      xlab="",
                      ylab="",
                      cexCol=0.1,
                      cexRow=0.4,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="target count",
                      cex.main=1)

  file <- paste0("../output/ tc.trend.2 ",dataset,"_",sigset," chem target pfas.xlsx")
  write.xlsx(res,file,rowNames=T)

  ##################################################################################
  map <- unique(mat[,c("dtxsid","name","category")])
  map <- map[is.element(map$name,rownames(res)),]
  map <- map[order(map$category),]
  res <- res[map$name,]
  names <- rownames(res)
  for(i in 1:nrow(res)) {
    name <- rownames(res)[i]
    category <- map[is.element(map$name,name),"category"]
    names[i] <- category
  }

  result <- heatmap.2(as.matrix(res),
                      margins=c(10,10),
                      dendrogram="col",
                      scale="none",
                      main="",
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
                      labRow=names)
  if(to.file) dev.off()
  else browser()
}
