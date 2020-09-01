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
#'
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#'  u2os_toxcast_pfas_pe1_normal
#'
#' Error bars are exp(er)*qt(.025,4) = exp(er)*2.7765
#--------------------------------------------------------------------------------------
pfas.tc.trend.3 <- function(to.file=F,
                            hccut=0.9,
                            tccut=1,
                            dbmdcut=30,
                            do.err.filter=T,
                            dataset="heparg2d_toxcast_pfas_pe1_normal",
                            sigcatalog="signatureDB_master_catalog 2020-08-14",
                            sigset="screen_large",
                            method="fc",
                            do.load=F) {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/PFAS/pfas.trend.3 ",dataset,"_",sigset,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  if(do.load) {
    annotations <- signatureCatalogLoader(sigset,sigcatalog)
    ANN <<- annotations
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat <- SIGNATURE_CR
    MAT <<- mat
  }
  mat <- MAT
  mat[is.na(mat$super_target),"super_target"] <- "-"
  catalog <- ANN
  file <- "../input/pfas/QC sample map 2020-05-04.xlsx"
  qc <- read.xlsx(file)
  #qc <- qc[is.element(qc$score,c("H","M","No score","TBD")),]
  spid.list <- qc$spid
  mat <- mat[is.element(mat$sample_id,spid.list),]

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

  vals <- c(1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6)
  nval <- length(vals)
  x <- NULL
  y <- NULL
  for(val in vals) {
    bottom <- val
    top <- val+0.5
    temp <- mat[mat$top_over_cutoff>=bottom,]
    temp <- temp[temp$top_over_cutoff<top,]
    yy <- temp$hitcall
    xx <- vector(mode="character",length=length(yy))
    xx[] <- paste0(bottom,"-",top)
    y <- c(y,yy)
    x <- c(x,xx)
  }
  boxplot(y~x,cex.lab=1.2,cex.axis=1.2,xlab="top / cutoff",ylab="hitcall")

  vals <- c(0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95)
  nval <- length(vals)
  x <- NULL
  y <- NULL
  for(val in vals) {
    bottom <- val
    top <- val+0.05
    temp <- mat[mat$hitcall>=bottom,]
    temp <- temp[temp$hitcall<top,]
    yy <- temp$bmd
    xx <- vector(mode="character",length=length(yy))
    xx[] <- paste0(bottom,"-",top)
    y <- c(y,yy)
    x <- c(x,xx)
  }
  boxplot(y~x,log="y",cex.lab=1.2,cex.axis=1.2,xlab="hitcall",ylab="BMD")


  ##################################################################################
  hc <- mat$hitcall
  tc <- mat$top_over_cutoff
  bmd <- mat$bmd
  er <- mat$er
  ci <- exp(er)*2.7765
  delta.bmd <- mat$bmdu-mat$bmdl
  delta.bmd <- delta.bmd[!is.na(delta.bmd)]
  hist(hc)
  hist(tc)
  hist(ci)

  err.range <- exp(er)*2.7765
  err_over_cutoff <- err.range / mat$cutoff
  hist(err_over_cutoff)

  #lbmd <- log10(bmd)
  #model <- lm(lbmd~tc+hc)
  #sm <- summary(model)
  #coefs <- sm$coefficients
  #intercept <- coefs[1,1]
  #ctc <- coefs[2,1]
  #chc <- coefs[3,1]
  #lbmdpred <- intercept + ctc*tc + chc*hc
  #plot(lbmdpred~lbmd,pch=".")

  #plot(tc~err_over_cutoff,pch=".",main="all",xlim=c(0,8),ylim=c(0,8))
  #lines(c(0,8),c(0,8))
  #y <- tc[hc>0.9]
  #x <- err_over_cutoff[hc>0.9]
  #plot(y~x,pch=".",main="hc>0.9",xlab="err_over_cutoff",ylab="tc",xlim=c(0,8),ylim=c(0,8))
  #lines(c(0,8),c(0,8))

  ##################################################################################
  # Do the filtering
  do.filter <- F
  if(do.filter) {
    mat <- MAT
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
    cat("allcut:",format(nrow(mat)/nrow(MAT),digits=3),"\n")

    cat("=============================================================\n")

    mat <- mat[is.element(mat$sample_id,spid.list),]
  }

  ##################################################################################
  # develop the score matrix
  mat <- MAT
  mat <- mat[is.element(mat$sample_id,spid.list),]
  mat <- mat[mat$hitcall>hccut,]
  mat <- mat[mat$top_over_cutoff>tccut,]
  ##################################################################################
  par(mfrow=c(3,2),mar=c(5,6,4,4))
  # Look at the simple trends
  file <- "../input/pfas/PFAS info.xlsx"
  info <- read.xlsx(file)
  info <- info[!is.na(info$endgroup),]

  for(nc in c(6,7,8)) {
    temp <- info[info$nc==nc,]
    g.list <- unique(temp$endgroup)
    x <- NULL
    y <- NULL
    for(g in g.list) {
      dtxsid.list <- unique(temp[is.element(temp$endgroup,g),"dtxsid"])
      bvals <- mat[is.element(mat$dtxsid,dtxsid.list),"bmd"]
      xvals <- vector(mode="character",length=length(bvals))
      xvals[] <- g
      x <- c(x,xvals)
      y <- c(y,bvals)
    }
    boxplot(y~x,main=paste("CF: ",nc),ylab="",xlab="BMD",log="x",ylim=c(0.001,100),
            horizontal=T, las=1,cex.lab=0.8,cex.axis=0.8)
  }

  for(g in c("telomer alcohol","sulfonic acid","alcohol","amide","acid","sulfonate","phosphonic acid")) {
    temp <- info[is.element(info$endgroup,g),]
    nc.list <- unique(temp$nc)
    x <- NULL
    y <- NULL
    for(nc in nc.list) {
      dtxsid.list <- unique(temp[temp$nc==nc,"dtxsid"])
      bvals <- mat[is.element(mat$dtxsid,dtxsid.list),"bmd"]
      xvals <- vector(mode="character",length=length(bvals))
      if(g=="acid") {
        cnc <- as.character(nc)
        if(nc<10) cnc <- paste0("0",cnc)
        nc <- cnc
      }
      xvals[] <- nc
      x <- c(x,xvals)
      y <- c(y,bvals)
    }
    cat(g,length(x),"\n")
    if(length(x)>0)
      boxplot(y~x,main=paste("End group: ",g),ylab="",xlab="BMD",log="x",ylim=c(0.001,100),
              horizontal=T, las=1)
  }

  ##################################################################################

  st.list <- sort(unique(ANN$super_target))

  st.list <- st.list[!is.element(st.list,c("-","unknown"))]
  res1 <- unique(mat[,c("sample_id","dtxsid","casrn","name")])
  res1 <- res1[order(res1$dtxsid),]
  rownames(res1) <- res1$sample_id
  dtxsid.list <- res1$dtxsid
  dups <- duplicated(dtxsid.list)
  temp <- dtxsid.list[dups]
  res1$duplicate <- 0
  res1[is.element(res1$dtxsid,temp),"duplicate"] <- 1

  nchem <- nrow(res1)
  res2 <- matrix(nrow=nchem,ncol=length(st.list))
  colnames(res2) <- st.list
  rownames(res2) <- rownames(res1)
  res2[] <- 0
  for(st in st.list) {
    temp <- ANN[is.element(ANN$super_target,st),]
    denom <- nrow(temp)
    temp <- mat[is.element(mat$super_target,st),]
    if(nrow(temp)>0) {
      cat(st,":",nrow(temp),"\n")
      for(i in 1:nrow(temp)) {
        sid <- temp[i,"sample_id"]
        bmd <- temp[i,"bmd"]
        lbmd <- -log10(bmd/1000)
        score <- abs(temp[i,"top"]) * lbmd / denom
        res2[sid,st] <- res2[sid,st] + score
      }
    }
  }
  cs <- colSums(res2)
  res2 <- res2[,cs>1.5]
  hist(cs)

  res3 <- cbind(res1,res2)
  file <- paste0("../output/PFAS/pfas_httr_score ",dataset,"_",sigset,".xlsx")
  write.xlsx(res3,file)

  colors <- res1$name
  rn <- colors
  dups <- duplicated(rn)
  dups <- rn[dups]
  colors[] <- "black"
  colors[is.element(rn,dups)] <- "red"

  result <- heatmap.2(as.matrix(res2),
                      margins=c(10,10),
                      dendrogram="both",
                      scale="none",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.4,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="target count",
                      cex.main=1,
                      labRow=res1[,"name"],
                      main="All PFAS",
                      RowSideColors=colors,
                      colRow=colors)

  res2a <- res2[res1$duplicate==1,]
  res1a <- res1[res1$duplicate==1,]
  result <- heatmap.2(as.matrix(res2a),
                      margins=c(20,20),
                      dendrogram="col",
                      scale="none",
                      xlab="",
                      ylab="",
                      cexCol=0.1,
                      cexRow=1,
                      Rowv=F,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="target count",
                      cex.main=1,
                      labRow=res1a[,"name"],
                      main="Duplicated PFAS")
  ############################################################################
  # heat map by category
  file <- "../input/pfas/PFAS categories.xlsx"
  cats <- read.xlsx(file)
  cat.list <- sort(unique(cats$category))
  for(categ in cat.list) {
    dtxsid.list <- cats[is.element(cats$category,categ),"dtxsid"]
    if(length(dtxsid.list)>2) {
      res2a <- res2[is.element(res1$dtxsid,dtxsid.list),]
      res1a <- res1[is.element(res1$dtxsid,dtxsid.list),]
      if(nrow(res1a)>2) {
        cat(categ,length(dtxsid.list),nrow(res1a),nrow(res2a),"\n")
        result <- heatmap.2(as.matrix(res2a),
                            margins=c(10,20),
                            dendrogram="both",
                            scale="none",
                            xlab="",
                            ylab="",
                            cexCol=0.1,
                            cexRow=0.8,
                            Rowv=T,
                            Colv=T,
                            trace="none",
                            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                            key=T,
                            col=brewer.pal(9,"Reds"),
                            key.title="Key",
                            key.xlab="target count",
                            cex.main=1,
                            labRow=res1a[,"name"],
                            main=categ,
                            cex.main=0.5)
      }
    }
  }
  #  box plot of scores by category
  #  score sume ~ category, mw
  # most common signatures or super targets
  # model score ~ category, mw


  if(to.file) dev.off()
  else browser()

}
