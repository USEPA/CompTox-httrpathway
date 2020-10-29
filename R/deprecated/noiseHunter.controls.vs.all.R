#--------------------------------------------------------------------------------------
#' Make heatmaps of the FCMAT and FCMAT.SE by plate groups
#'
#' @param to.file If TRUE, send plots to a pdf
#' @param do.load If TRUE, load the data into a global
#' @param dataset The name of the data set to analyze
#' @param sigset The signature set to use
#' @param sigcatalog The signature catalog to use
#' @param method The signature scoring method
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
noiseHunter.controls.vs.all <- function(to.file=F,
                                        do.load=F,
                                        dataset="DMEM_6hr_screen_normal_pe_1",
                                        dir="../input/httr_mcf7_screen/meanncnt0_5-plateteffect_0-shrinkage_normal_DMEM_6_controls/"){
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter.controls.pdf")
     pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,4),mar=c(4,3,3,2))

  if(do.load) {
    file <- paste0("../output/noiseHunter/noiseHunter1 ",dataset," pg distribution.xlsx")
    mvals <- read.xlsx(file)
    MVALS <<- mvals
  }
  mvals <- MVALS
  prefix <- "httr_mcf7_ph1_"
  type.list <- c("QC_BL_TSA","QC_HBRR-vs-UHRR","refChem_GEN","refChem_SIRO","refChem_TSA")
  type.list <- c("QC_BL_TSA","refChem_TSA")
  type.list <- c("QC_BL_TSA","refChem_GEN","refChem_SIRO","refChem_TSA")
  suffix <- "_meanncnt0_5-plateteffect_1-shrinkage_normal_fc.tsv"
  for(pg in 1:48) {
    mat <- NULL
    for(type in type.list) {
      file <- paste0(dir,prefix,type,"_pg",pg,suffix)
      print(type)
      temp <- read.table(file,sep="\t",stringsAsFactors=F,header=T,quote="")
      print(names(temp))
      if(contains(type,"refChem")) {
        name.list <- names(temp)
        name.list <- name.list[!is.element(name.list,"chem_id")]
        temp <- temp[,name.list]
      }
      mat <- rbind(mat,temp)
    }
    print(dim(mat))
    probes <- mat$probe_id
    genes <- separate(mat,"probe_id",sep="_",into=c("genes","id"))
    mat$gene <- genes[,1]
    print(dim(mat))

    res <- reshape2::dcast(mat,gene~trt_name,value.var="log2FoldChange",fill=0,fun.aggregate=max)
    rownames(res) <- res[,1]
    res <- as.matrix(res[,2:ncol(res)])
    mval.pg <- mvals[mvals$pg==pg,]
    rownames(mval.pg) <- mval.pg$gene
    gene.list <- mval.pg$gene
    gene.list <- gene.list[is.element(gene.list,rownames(res))]
    res <- res[gene.list,]
    mval.pg <- mval.pg[gene.list,]
    mval.pg <- mval.pg[order(mval.pg$mean),]
    gene.list <- mval.pg$gene
    res <- res[gene.list,]

    x <- mval.pg$mean
    do.plot1 <- F
    if(do.plot1) {
      y <- mval.pg$median
      xlim <- 0.4
      ylim <- 4
      plot(y~x,cex.lab=1.2,cex.axis=1.2,main=paste("PG:",pg),
           xlim=c(-xlim,xlim),ylim=c(-ylim,ylim),
           xlab="mean(l2fc)",ylab="l2fc(comparator)",type="l")
      col.list <- c("red","black")
      xval <- -xlim
      yval <- -2
      for(i in 1:ncol(res)) {
        y <- res[,i]
        points(y~x,pch=".",col=col.list)
        result <- lm(y~x)
        #print(summary(result))
        pval <- summary(result)$coefficients[2,4]
        text(xval,yval,format(pval,digits=2),pos=4)
        yval <- yval-1
        #browser()
      }
    }
    ylims <- c(3,3,3,3)
    for(i in 1:ncol(res)) {
      y <- res[,i]
      sp <- cor.test(x,y,method="spearman")
      rho <- sp$estimate
      paste("PG:",pg,"\n",type.list[i])
      main <- paste(type.list[i],"\nPG:",pg," rho=",format(rho,digits=2))
      nrow <- length(y)
      index <- y
      index[] <- 1
      index[(nrow/4):(nrow/2)] <- 2
      index[(nrow/2):(3*nrow/4)] <- 3
      index[(3*nrow/4):(nrow)] <- 4
      ylim <- ylims[i]
      boxplot(y~index,main=main,
              cex.main=1,cex.lab=1.2,cex.axis=1.2,ylim=c(-ylim,ylim),
              xlab="|l2fc| quartile",ylab="control l2fc")
      lines(c(-100,100),c(0,0))
    }


    nrow <-

    if(!to.file) browser()
  }
  if(to.file) dev.off()
}
