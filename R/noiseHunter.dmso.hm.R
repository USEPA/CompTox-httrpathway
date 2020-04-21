#--------------------------------------------------------------------------------------
#' Find issues with the DMSO samples that may be causing the odd noise issue
#'
#' @param dir The directory to work in
#'
#--------------------------------------------------------------------------------------
noiseHunter.dmso.hm <- function(to.file=F,
                                do.read=F,
                                dirin="../input/httr_mcf7_screen/raw_ctrl/",
                                welltype="dmso",
                                well_type_prefix="",
                                dirout="../output/noiseHunter/",
                                do.prep=F,
                                do.1=F,
                                do.2=F,
                                do.3=F,
                                do.4=F,
                                cutoff=5,
                                ngene=-1) {
  printCurrentFunction()
  if(do.read) {
    file <- paste0(dirin,"httr-ph-i-raw-pl-bl_",welltype,"-v1.tsv")
    print(file)
    mat <- read.delim(file,stringsAsFactors=F)
    x <- unite(mat[,c("plate_id","block_id","pg_id","well_id")],"rowkey")
    mat$rowkey <- x[,1]
    CTRLMAT <<- mat
    pmat <- unique(CTRLMAT[,c("plate_id","block_id","pg_id","well_id")])
    x <- unite(pmat,"rowkey")
    rownames(pmat) <- x[,1]
    PLATEMAT <<- pmat
  }

  if(do.1) {
    if(to.file) {
      fname <- paste0(dirout,"DMSO raw counts.pdf")
      pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(5,4),mar=c(4,4,2,2))

    pg.list <- sort(unique(CTRLMAT$pg_id))
    for(pg in pg.list) {
      mat <- CTRLMAT[is.element(CTRLMAT$pg_id,pg),]
      x <- mat$probe_count
      lx <- log2(x)
      hist(lx,main=paste("PG:",pg),xlab="log2(raw counts)",
           cex.axis=1.2,cex.lab=1.2)
      mf <- mat$mapd_frac
      breaks <- seq(from=0.6,to=1.,by=0.01)
      hist(mf,main=paste("PG:",pg),xlab="mapped fraction",
           cex.axis=1.2,cex.lab=1.2,breaks=breaks)
      if(!to.file) browser()
    }
    if(to.file) dev.off()
  }
  if(do.2) {
    if(to.file) {
      fname <- paste0(dirout,"DMSO principal components 2.pdf")
      pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(2,1),mar=c(4,4,2,2))
    pg.list <- sort(unique(CTRLMAT$pg_id))
    if(!exists("RES")) {
      newmat <- NULL
      for(pg in pg.list) {
        mat <- CTRLMAT[is.element(CTRLMAT$pg_id,pg),]
        sid.list <- unique(mat$sample_id)
        mat$spidpgid <- "-"
        for(i in 1:length(sid.list) ) {
          sample_id <- sid.list[i]
          spidpgid <- paste0("pg_",pg,"_",i)
          mat[mat$sample_id==sample_id,"spidpgid"] <- spidpgid
        }
        newmat <- rbind(newmat,mat)
      }
      res <- reshape2::dcast(newmat,spidpgid~probe_id,value.var="probe_count",fill=0.5)
      rownames(res) <- res[,1]
      res <- res[,2:ncol(res)]
      res <- as.matrix(res)
      rs <- rowSums(res)
      for(i in 1:nrow(res)) res[i,] <- res[i,]*1000000/rs[i]
      res <- log2(res)
      RES <<- res
    }
    res <- RES
    pres <- prcomp(res)
    pc <- pres$x
    col.list.0 <- c("red","orange","yellow","green",
                    "blue","violet","gray","black",
                    "white","cyan","pink","khaki")
    pch.list.0 <- c(21,22,23,24)
    col.list <- vector(length=length(pg.list),mode="character")
    col.list[] <- "white"
    pch.list <- vector(length=length(pg.list),mode="integer")
    pch.list[] <- 1
    pg <- 1
    for(i in 1:12) {
      for(j in 1:4) {
        col.list[pg] <- col.list.0[i]
        pch.list[pg] <- pch.list.0[j]
        pg <- pg+1
      }
    }
    plot(pc[,"PC2"]~pc[,"PC1"],xlim=c(-200,200),ylim=c(-150,150),
         main=paste("MCF7 DMSO"),cex.lab=1.2,cex.axis=1.2,xlab="PC1",ylab="PC2",type="n")
    for(i in 1:nrow(pc)) {
      x <- pc[i,"PC1"]
      y <- pc[i,"PC2"]
      rn <- rownames(pc)[i]
      rn <- substr(rn,4,5)
      pg <- as.numeric(str_replace(rn,"_",""))
      color <- col.list[pg]
      pch <- pch.list[pg]
      points(x,y,pch=pch,bg=color)
      if(pg==7) text(x,y,paste0(pg),pos=4)
    }

    xy <- pc[,c("PC1","PC2","PC3","PC4")]
    colnames(xy)[4] <- "pg"
    for(i in 1:nrow(xy)) {
      rn <- rownames(xy)[i]
      rn <- substr(rn,4,5)
      pg <- as.numeric(str_replace(rn,"_",""))
      xy[i,"pg"] <- pg
    }

    for(i in 1:48) {
      temp <- xy[xy[,"pg"]==i,c("PC1","PC2")]
      x <- dist(temp)
      cat("pg: ",i, max(as.numeric(x)),"\n")
    }
    if(!to.file) browser()

    plot(pc[,"PC3"]~pc[,"PC2"],xlim=c(-100,100),ylim=c(-100,100),
         main=paste("MCF7 DMSO"),cex.lab=1.2,cex.axis=1.2,xlab="PC2",ylab="PC3",type="n")
    for(i in 1:nrow(pc)) {
      x <- pc[i,"PC2"]
      y <- pc[i,"PC3"]
      rn <- rownames(pc)[i]
      rn <- substr(rn,4,5)
      pg <- as.numeric(str_replace(rn,"_",""))
      color <- col.list[pg]
      pch <- pch.list[pg]
      points(x,y,pch=pch,bg=color)
      if(pg==7) text(x,y,paste0(pg),pos=4)
    }


    if(!to.file) browser()

    if(to.file) dev.off()
  }

  if(do.3) {
    if(to.file) {
      fname <- paste0(dirout,"DMSO principal components HM.pdf")
      pdf(file=fname,width=7,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    pg.list <- sort(unique(CTRLMAT$pg_id))
    pg.list <- c(4,7)
    if(!exists("RES")) {
      newmat <- NULL
      for(pg in pg.list) {
        mat <- CTRLMAT[is.element(CTRLMAT$pg_id,pg),]
        sid.list <- unique(mat$sample_id)
        mat$spidpgid <- "-"
        for(i in 1:length(sid.list) ) {
          sample_id <- sid.list[i]
          spidpgid <- paste0("pg_",pg,"_",i)
          mat[mat$sample_id==sample_id,"spidpgid"] <- spidpgid
        }
        newmat <- rbind(newmat,mat)
      }
      res <- reshape2::dcast(newmat,spidpgid~probe_id,value.var="probe_count",fill=0.5)
      rownames(res) <- res[,1]
      res <- res[,2:ncol(res)]
      res <- as.matrix(res)
      rs <- rowSums(res)
      for(i in 1:nrow(res)) res[i,] <- res[i,]*1000000/rs[i]
      res <- log2(res)
      RES <<- res
    }
    res <- RES
    if(ngene>0) res <- res[,1:ngene]
    result <- heatmap.2(res,
                        margins=c(5,5),
                        dendrogram="both",
                        scale="none",
                        main=paste("DMSO CPM"),
                        xlab="",
                        ylab="",
                        cexCol=0.1,
                        cexRow=1,
                        Rowv=T,
                        Colv=T,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(9,"Reds"),
                        key.title="Key",
                        key.xlab="CPM",
                        cex.main=1)



    if(!to.file) browser()
    else dev.off()
  }

}
