#source("https://bioconductor.org/biocLite.R")
#.libPaths( "C:/Program Files/R/R-3.3.3/library")
library(fastmatch)
library(reshape2)
#source("R/general_utils.R")
#options(java.parameters = "-Xmx8000m")
#library(httRws)
library(openxlsx)
library(tidyr)
library(grDevices)
library(RColorBrewer)
library(gplots)
library(ape)
library(rawr)
#--------------------------------------------------------------------------------------
#' Prepare the raw data for the replicated chemicals
#'
#' @param dir The directory to work in
#'
#--------------------------------------------------------------------------------------
raw.data.analysis <- function(to.file=F,
                              welltype="bl_tsa",
                              do.read.dmso=F,
                              do.read.other=F,
                              prep.l2fcmat=F,
                              analysis_id=3,
                              read_id=-1,
                              welltypeA="hbrr",
                              welltypeB="uhrr",
                              cutoff=0.95,
                              l2fcmax=10,
                              dirin="analysis/raw_count_analysis/data/",
                              dirout="analysis/raw_count_analysis/") {
  printCurrentFunction()
  ############################################################################################
  ############################################################################################
  ############################################################################################
  if(do.read.dmso) {
    cat("read DMSO\n")
    wt <- "vehicle_control"
    file <- paste0(dirin,"data_",wt,"_",cutoff,".RData")
    print(file)
    load(file)
    PLATEMAT.DMSO <<- PLATEMAT
    COUNTMAT.DMSO <<- COUNTMAT
  }
  ############################################################################################
  ############################################################################################
  ############################################################################################
  if(do.read.other) {
    cat("read ",welltype,"\n")
    file <- paste0(dirin,"data_",welltype,"_",cutoff,".RData")
    print(file)
    load(file)
    PLATEMAT <<- PLATEMAT
    COUNTMAT <<- COUNTMAT
  }
  ############################################################################################
  ############################################################################################
  ############################################################################################
  if(prep.l2fcmat) {
    cat("create l2fcmat ",welltype,"\n")
    countmat.dmso <- COUNTMAT.DMSO
    countmat <- COUNTMAT
    probe.list <- sort(colnames(countmat.dmso))
    probe.list <- probe.list[is.element(probe.list,colnames(countmat))]
    countmat.dmso <- countmat.dmso[,probe.list]
    countmat <- countmat[,probe.list]
    plate.list <- sort(unique(PLATEMAT.DMSO[,"plate_id"]))
    l2fcmat <- matrix(nrow=nrow(countmat),ncol=length(probe.list))
    colnames(l2fcmat) <- probe.list
    rownames(l2fcmat) <- rownames(countmat)
    for(i in 1:length(plate.list)) {
      plate <- plate.list[i]
      sample.list.dmso <- rownames(PLATEMAT.DMSO[is.element(PLATEMAT.DMSO$plate_id,plate),])
      meds.dmso <- apply(countmat.dmso[sample.list.dmso,],FUN=median,MARGIN=2)
      sample.list <- rownames(PLATEMAT[is.element(PLATEMAT$plate_id,plate),])
      for(sample in sample.list) {
        scount <- countmat[sample,]
        l2fcmat[sample,] <- scount-meds.dmso
      }
    }
    file <- paste0(dirin,"l2fcmat_",welltype,"_",cutoff,".RData")
    L2FCMAT <- l2fcmat
    save(L2FCMAT,file=file)
  }

  ############################################################################################
  ############################################################################################
  ############################################################################################
  if(analysis_id==1) {
    cat("run analysis 1\n")
    if(to.file) {
      fname <- paste0(dirout,"raw_data_analysis_dmso_correlation.pdf")
      pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(3,2),mar=c(4,4,2,2))
    countmat <- COUNTMAT.DMSO
    meds <- apply(countmat,FUN=median,MARGIN=2)
    platemat <- PLATEMAT.DMSO
    name.list <- c("sample_id","plate_id","block_id","pg_id","well_id","mad")
    madmat <- as.data.frame(matrix(nrow=nrow(platemat),ncol=length(name.list)))
    names(madmat) <- name.list
    madmat[,"plate_id"] <- platemat[,"plate_id"]
    madmat[,"pg_id"] <- platemat[,"pg_id"]
    madmat[,"block_id"] <- platemat[,"block_id"]
    madmat[,"well_id"] <- platemat[,"well_id"]
    for(i in 1:nrow(platemat)) {
      sample_id <- rownames(platemat)[i]
      madmat[i,"sample_id"] <- sample_id
      cat(sample_id,"\n")
      vals <- countmat[sample_id,]
      delta <- vals-meds
      mval <- mad(delta)
      madmat[i,"mad"] <- mval
      plot(vals~meds,main=sample_id,xlab="log2(DMSO median)",ylab="log2(DMSO well)",
           xlim=c(0,15),ylim=c(0,15),pch=21,cex=0.2)
      lines(c(0,100),c(0,100))
      text(0,14,paste("mad: ",format(mval,digits=3)),pos=4)
      if(!to.file) browser()
    }
    if(to.file) dev.off()
    file <- paste0(dirout,"raw_data_analysis_dmso_correlation_mad.xlsx")
    write.xlsx(madmat,file)
  }
  ############################################################################################
  ############################################################################################
  ############################################################################################
  if(analysis_id==2) {
    cat("run analysis 2 ",welltype,"\n")

    if(to.file) {
      fname <- paste0(dirout,"raw_data_analysis_l2fc_correlation_",welltype,".pdf")
      pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(3,2),mar=c(4,4,2,2))

    file <- paste0(dirin,"l2fcmat_",welltype,"_",cutoff,".RData")
    load(file)
    l2fcmat <- L2FCMAT

    plate.list <- rownames(l2fcmat)
    meds <- apply(l2fcmat,FUN=median,MARGIN=2)
    sample.list <- rownames(l2fcmat)
    x <- str_split(sample.list,"_")
    madmat <- as.data.frame(matrix(unlist(x),nrow=length(sample.list),byrow=T))
    madmat <- cbind(sample.list,madmat)
    names(madmat) <- c("sample_id","plate_id","block_id","pg_id","well_id")
    madmat$mad <- NA

    for(sample_id in sample.list) {
      cat(sample_id,"\n")

      delta <- l2fcmat[sample_id,]-meds
      mval <- mad(delta)
      madmat[is.element(madmat[,"sample_id"],sample_id),"mad"] <- mval

      plot(l2fcmat[sample_id,]~meds,main=sample_id,
           xlab=paste0("l2fc(",welltype,"|median)"),
           ylab=paste0("l2fc(",welltype,"|plate)"),
           pch=21,cex=0.2,xlim=c(-l2fcmax,l2fcmax),ylim=c(-l2fcmax,l2fcmax))
      lines(c(-100,100),c(-100,100))
      lines(c(-100,100),c(0,0))
      lines(c(0,0),c(-100,100))
      text(-l2fcmax,l2fcmax,paste("mad: ",format(mval,digits=3)),pos=4)
      if(!to.file) browser()
    }
    file <- paste0(dirout,"raw_data_analysis_l2fc_correlation_mad_",welltype,".xlsx")
    write.xlsx(madmat,file)

    if(to.file) dev.off()
  }

  ############################################################################################
  ############################################################################################
  ############################################################################################
  if(read_id==3) {
    cat("read for analysis 3 \n")
    file <- paste0(dirin,"data_",welltypeA,"_",cutoff,".RData")
    print(file)
    load(file)
    PLATEMAT.A <<- PLATEMAT
    COUNTMAT.A <<- COUNTMAT
    file <- paste0(dirin,"data_",welltypeB,"_",cutoff,".RData")
    print(file)
    load(file)
    PLATEMAT.B <<- PLATEMAT
    COUNTMAT.B <<- COUNTMAT
  }
  if(analysis_id==3) {
    cat("run analysis 3 \n")

    if(to.file) {
      fname <- paste0(dirout,"raw_data_analysis_l2fc_correlation_",welltypeA,"_",welltypeB,".pdf")
      pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(3,2),mar=c(4,4,2,2))

    matA <- as.matrix(COUNTMAT.A)
    matB <- as.matrix(COUNTMAT.B)
    gene.list <- colnames(matA)
    gene.list <- gene.list[is.element(gene.list,colnames(matB))]
    matA <- matA[,gene.list]
    matB <- matB[,gene.list]
    plateA <- PLATEMAT.A
    plateB <- PLATEMAT.B
    plate.list <- rownames(plateA)

    medA <- apply(matA,FUN=median,MARGIN=2)
    medB <- apply(matB,FUN=median,MARGIN=2)
    l2fcG <- medA-medB
    plate.list <- unique(plateA[,"plate_id"])

    madmat <- unique(plateA[,c("plate_id","block_id","pg_id")])
    madmat$mad <- NA
    rownames(madmat) <- madmat[,"plate_id"]

    for(plate in plate.list) {
      cat(plate,"\n")
      mA <- apply(matA[rownames(plateA[is.element(plateA[,"plate_id"],plate),]),],FUN=median,MARGIN=2)
      mB <- apply(matB[rownames(plateB[is.element(plateB[,"plate_id"],plate),]),],FUN=median,MARGIN=2)
      l2fc <- mA-mB

      delta <- l2fc-l2fcG
      mval <- mad(delta)
      madmat[plate,"mad"] <- mval

      plot(l2fc~l2fcG,main=plate,
           xlab=paste0("l2fc(",welltypeA,"-",welltypeB,"|global)"),
           ylab=paste0("l2fc(",welltypeA,"-",welltypeB,"|plate)"),
           pch=21,cex=0.2,xlim=c(-l2fcmax,l2fcmax),ylim=c(-l2fcmax,l2fcmax))
      lines(c(-100,100),c(-100,100))
      lines(c(-100,100),c(0,0))
      lines(c(0,0),c(-100,100))
      text(-l2fcmax,l2fcmax,paste("mad: ",format(mval,digits=3)),pos=4)
      if(!to.file) browser()
    }

    file <- paste0(dirout,"raw_data_analysis_l2fc_correlation_mad_",welltypeA,"_",welltypeB,".xlsx")
    write.xlsx(madmat,file)
    if(to.file) dev.off()
  }
  ############################################################################################
  ############################################################################################
  ############################################################################################
  if(analysis_id==4) {
    cat("run analysis 4 \n")
    welltype.list <- c("bl_dmso","bl_tsa","hbrr","uhrr","untreated","vehicle_control","TSA","GEN","SIRO")
    #welltype.list <- c("vehicle_control","TSA","GEN","SIRO")
    ngene <- 50000
    ngene <- 10000
    for(i in 1:length(welltype.list)) {
      welltype <- welltype.list[i]
      if(to.file) {
        fname <- paste0(dirout,"raw_data_analysis_heatmap_",welltype,".pdf")
        pdf(file=fname,width=7,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
      }
      file <- paste0(dirin,"data_",welltype,"_",cutoff,".RData")
      print(file)
      load(file)
      PLATEMAT <<- PLATEMAT
      COUNTMAT <<- COUNTMAT
      gene.list <- colnames(COUNTMAT)
      if(i==1) {
        if(ngene<length(gene.list)) gene.list.use <- sample(gene.list,ngene)
        else gene.list.use <- gene.list
      }
      gene.list <- gene.list.use[is.element(gene.list.use,gene.list)]
      #browser()
      mat <- COUNTMAT[,gene.list]
      block.list <- PLATEMAT[rownames(COUNTMAT),"block_id"]
      mat <- mat[block.list<5,]
      col.list <- PLATEMAT[rownames(mat),"block_id"]
      col.list[is.element(col.list,1)] <- "red"
      col.list[is.element(col.list,2)] <- "black"
      col.list[is.element(col.list,3)] <- "white"
      col.list[is.element(col.list,5)] <- "cyan"
      mat[mat<0] <- 0
      cs <- colMed(mat)
      #browser()
      for(k in 1:10) {
        mat <- rbind(cs,mat)
        rownames(mat)[1] <- "medium"
        col.list <- c("green",col.list)
      }
      cat(welltype,"start heatmap\n")
      heatmap.3(mat,
                Rowv=T,
                Colv=T,
                #hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                #dendrogram="both",
                symm=F,
                scale="none",
                cols=brewer.pal(9,"Reds"),
                margins=c(5,5),
                cexCol=0.1,
                cexRow=0.1,
                key=T,
                main=paste("log2(counts): ",welltype),
                xlab="",
                ylab="",
                trace="none",
                key.title="Key",
                #key.xlab="log2(counts)",
                cex.main=1,
                RowSideColors=col.list)
      if(!to.file) browser()
      else dev.off()
    }
  }
  ############################################################################################
  ############################################################################################
  ############################################################################################
  if(analysis_id==5) {
    cat("run analysis 5 \n")
    welltype.list <- c("bl_dmso","bl_tsa","hbrr","uhrr","untreated","vehicle_control","TSA","GEN","SIRO")
    #welltype.list <- c("vehicle_control","TSA","GEN","SIRO")
    #welltype.list <- c("vehicle_control","TSA")
    #welltype.list <- c("untreated","vehicle_control","TSA","GEN","SIRO")
    ngene <- 50000
    ngene <- 100
    mad.list <- vector("list", length(welltype.list))
    names(mad.list) <- welltype.list
    med.list <- mad.list

    gene.list.all <- NULL
    for(i in 1:length(welltype.list)) {
      welltype <- welltype.list[i]
      file <- paste0(dirin,"data_",welltype,"_",cutoff,".RData")
      print(file)
      load(file)
      PLATEMAT <<- PLATEMAT
      COUNTMAT <<- COUNTMAT
      gene.list <- colnames(COUNTMAT)
       gene.list.all <- c(gene.list.all,gene.list)
      mat <- COUNTMAT[,gene.list]
      block.list <- PLATEMAT[rownames(COUNTMAT),"block_id"]
      mat <- mat[block.list<5,]
      mat[mat<0] <- 0
      col.med <- colMed(mat)
      col.mad <- colMad(mat)
      mad.list[[welltype]] <- col.mad
      med.list[[welltype]] <- col.med
    }

    gene.list <- sort(unique(gene.list))
    mat.mad <- as.data.frame(matrix(nrow=length(welltype.list),ncol=length(gene.list)))
    names(mat.mad) <- gene.list
    rownames(mat.mad) <- welltype.list
    for(i in 1:length(welltype.list)) {
      welltype <- welltype.list[i]
      temp <- mad.list[[welltype]]
      for(j in 1:length(temp)) {
        gene <- names(temp)[j]
        mat.mad[welltype,gene] <- temp[j]
      }
    }
    file <- paste0(dirout,"raw_data_analysis_heatmap_mad.xlsx")
    write.xlsx(t(mat.mad),file,rowNames=T,colNames=T)

    mat.med <- as.data.frame(matrix(nrow=length(welltype.list),ncol=length(gene.list)))
    names(mat.med) <- gene.list
    rownames(mat.med) <- welltype.list
    for(i in 1:length(welltype.list)) {
      welltype <- welltype.list[i]
      temp <- med.list[[welltype]]
      for(j in 1:length(temp)) {
        gene <- names(temp)[j]
        mat.med[welltype,gene] <- temp[j]
      }
    }
    file <- paste0(dirout,"raw_data_analysis_heatmap_med.xlsx")
    write.xlsx(t(mat.med),file,rowNames=T,colNames=T)


    if(to.file) {
      fname <- paste0(dirout,"raw_data_analysis_heatmap_mad_med.pdf")
      pdf(file=fname,width=7,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    cat("start heatmap\n")
    mat.mad[is.na(mat.mad)] <- 0
    heatmap.3(mat.mad,
              Rowv=T,
              Colv=T,
              symm=F,
              scale="none",
              cols=brewer.pal(9,"Reds"),
              margins=c(5,5),
              cexCol=0.1,
              cexRow=1,
              key=T,
              main=paste("log2(counts): MAD"),
              xlab="",
              ylab="",
              trace="none",
              key.title="Key",
              cex.main=1)
    if(!to.file) browser()

    mat.med[is.na(mat.med)] <- 0
    heatmap.3(mat.med,
              Rowv=T,
              Colv=T,
               symm=F,
              scale="none",
              cols=brewer.pal(9,"Reds"),
              margins=c(5,5),
              cexCol=0.1,
              cexRow=1,
              key=T,
              main=paste("log2(counts): Median"),
              xlab="",
              ylab="",
              trace="none",
              key.title="Key",
              cex.main=1)
    if(!to.file) browser()
    else dev.off()
  }
  ############################################################################################
  ############################################################################################
  ############################################################################################
  if(analysis_id==6) {
    cat("run analysis 6 \n")
    file <- paste0(dirout,"raw_data_analysis_heatmap_mad.xlsx")
    mad <- read.xlsx(file)
    file <- paste0(dirout,"raw_data_analysis_heatmap_med.xlsx")
    med <- read.xlsx(file)
    welltype.list <- c("bl_dmso","bl_tsa","hbrr","uhrr","untreated","vehicle_control","TSA","GEN","SIRO")
    rownames(mad) <- mad[,1]
    rownames(med) <- med[,1]

    gene.list <- rownames(mad)
    gene.list <- gene.list[is.element(gene.list,rownames(med))]
    mad <- mad[gene.list,]
    med <- med[gene.list,]
    if(to.file) {
      fname <- paste0(dirout,"raw_data_analysis_scatterplot_mad_med.pdf")
      pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(3,2),mar=c(4,4,2,2))

    for(welltype in welltype.list) {
      x <- med[,welltype]
      y <- mad[,welltype]
      plot(y~x,xlab="median",ylab="MAD",main=welltype,pch=".")
      if(!to.file) browser()
    }
    if(to.file) dev.off()
  }
  ############################################################################################
  ############################################################################################
  ############################################################################################
  if(analysis_id==7) {
    cat("run analysis 7 \n")
    file <- paste0(dirout,"raw_data_analysis_heatmap_mad.xlsx")
    mad <- read.xlsx(file)
    rownames(mad) <- mad[,1]
    file <- paste0(dirout,"raw_data_analysis_heatmap_med.xlsx")
    med <- read.xlsx(file)
    rownames(med) <- med[,1]
    file <- paste0(dirout,"BSP_Whole_Transcriptome_v1_v2_compare.xlsx")
    pdiff <- read.xlsx(file)
    pdiff <- pdiff[pdiff$v1==1,]
    drop <- pdiff[pdiff$v2==0,"Probe_ID"]
    keep <- pdiff[pdiff$v2==1,"Probe_ID"]
    welltype.list <- c("bl_dmso","bl_tsa","hbrr","uhrr","untreated","vehicle_control","TSA","GEN","SIRO")

    # box plot of mad across well types, drop, keep
    if(to.file) {
      fname <- paste0(dirout,"raw_data_analysis_v1_v2_probe.pdf")
      pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(2,2),mar=c(4,8,2,2))

    class.list <- NULL
    value.list <- NULL
    for(kd in c("keep","drop")) {
      if(kd=="keep") probe.list <- keep
      if(kd=="drop") probe.list <- drop
      for(welltype in welltype.list) {
        temp1 <- mad[probe.list,welltype]
        temp2 <- vector(length=length(temp1),mode="character")
        temp2[] <- paste(welltype,kd)
        class.list <- c(class.list,temp2)
        value.list <- c(value.list,temp1)
      }
    }
    boxplot(value.list~class.list,main="MAD",ylab="",horizontal=T,las=1,cex.lab=0.8,cex.axis=0.8)


    class.list <- NULL
    value.list <- NULL
    for(kd in c("keep","drop")) {
      if(kd=="keep") probe.list <- keep
      if(kd=="drop") probe.list <- drop
      for(welltype in welltype.list) {
        temp1 <- med[probe.list,welltype]
        temp2 <- vector(length=length(temp1),mode="character")
        temp2[] <- paste(welltype,kd)
        class.list <- c(class.list,temp2)
        value.list <- c(value.list,temp1)
      }
    }
    boxplot(value.list~class.list,main="Median",ylab="",horizontal=T,las=1,cex.lab=0.8,cex.axis=0.8)

        if(!to.file) browser()
    if(to.file) dev.off()
  }

}
