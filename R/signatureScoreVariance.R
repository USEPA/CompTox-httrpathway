
library(e1071)
#--------------------------------------------------------------------------------------
#' Visualize the variance by signature - this supplants signatureVariance.R
#'
#' @param min.ngene Signatures will only be saved if the number of genes is >= this value
#' @param max.ngene Signatures will only be saved if the number of genes is <= this value
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
signatureScoreVariance <- function(to.file=F,
                                   do.load=F,
                                   sigset="pilot_large_all_CMAP",
                                   sigcatalog="signatureDB_master_catalog 2020-01-31",
                                   basedir="../input/fcdata/",
                                   dataset="DMEM_6hr_screen_normal_pe_1",
                                   method="mygsea"){
  printCurrentFunction()

  if(do.load) {
    file <- paste0(basedir,"FCMAT2_",dataset,".RData")
    print(file)
    load(file)
    file <- paste0(basedir,"CHEM_DICT_",dataset,".RData",sep="")
    print(file)
    load(file)
    rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]

    FCMAT2 <<- FCMAT2
    CHEM_DICT <<- CHEM_DICT
    cat("   load signature data\n")
    file <- paste0("../input/signatures/signatureDB_genelists.RData")
    cat("   ",file,"\n")
    load(file) #genelists
    genelists <<- genelists

    annotations <<- signatureCatalogLoader(sigset,sigcatalog)

    file <- paste("../output/signature_score_summary/signaturescoremat_",sigset,"_",dataset,"_",method,".RData",sep="")
    print(file)
    load(file=file)
    signaturescoremat <<- signaturescoremat

    file <- paste("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData",sep="")
    print(file)
    load(file=file)
    signaturecr <<- SIGNATURE_CR

  }

  file <- "../input/signatures/httr_mcf7_ph1_DMSO_by_gene.xlsx"
  gcounts <- read.xlsx(file)
  if(to.file) {
    fname <- paste0("../output/signature_corr/signatureScoreVariance.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,4),mar=c(4,4,2,2))

  parent.list <- unique(annotations$parent)

  parent.list <- parent.list[is.element(parent.list,signaturecr$signature)]
  #parent.list <- parent.list[1:10]
  np <- length(parent.list)
  cat("number of parents; ",np,"\n")
  FCMAT2[is.nan(FCMAT2)] <- 0

  name.list <- c("signature","ngene",
                 "l2cpm.min","l2cpm.max",
                 "l2cpm.mean","l2cpm.sd",
                 "l2cpm.med","l2cpm.mad",
                 "gene.mean","gene.sd",
                 "gene.skewness","gene.kurtosis",
                 "gene.min","gene.max",
                 "gene.med","gene.mad",
                 "score.mean","score.sd",
                 "score.skewness","score.kurtosis",
                 "score.min","score.max",
                 "score.med","score.mad",

                 "l2fccpm.mean","l2fccpm.sd",
                 "l2fccpm.skewness","l2fccpm.kurtosis",
                 "l2fccpm.min","l2fccpm.max",
                 "l2fccpm.med","l2fccpm.mad",

                 "nhit","nhit.lt.1","nhit.lt.0.1"

  )
  res <- as.data.frame(matrix(nrow=length(np),ncol=length(name.list)))
  names(res) <- name.list

  for(i in 1:np) {
    parent <- parent.list[i]
    cat(i,parent,"\n")
    flush.console()
    sigs <- annotations[is.element(annotations$parent,parent),"signature"]
    gl <- NULL
    for(sig in sigs) gl <- c(gl,genelists[sig][[1]])
    gl <- gl[is.element(gl,colnames(FCMAT2))]

    temp1 <- gcounts[is.element(gcounts$gene_symbol,gl),c("gene_symbol","log2_CPM_mean")]
    x <- temp1$log2_CPM_mean
    res[i,"l2cpm.mean"] <- mean(x)
    res[i,"l2cpm.sd"] <- sd(x)
    res[i,"l2cpm.med"] <- median(x)
    res[i,"l2cpm.mad"] <- mad(x)
    res[i,"l2cpm.min"] <- min(x)
    res[i,"l2cpm.max"] <- max(x)
    if(length(x)>5) plot(density(x),main=paste("l2cpm\n",parent),cex.lab=1.2,cex.axis=1.2,xlab="l2CPM",xlim=c(-1,15),cex.main=0.8)
    else plot(c(1,1),main="missing plot")

    temp1 <- temp1[!duplicated(temp1$gene_symbol),]
    rownames(temp1) <- temp1$gene_symbol

    temp2 <- FCMAT2[,gl]

    x <- as.numeric(temp2)
    x <- sort(x)
    res[i,"signature"] <- sig
    res[i,"ngene"] <- length(gl)
    res[i,"gene.mean"] <- mean(x)
    res[i,"gene.sd"] <- sd(x)
    res[i,"gene.med"] <- median(x)
    res[i,"gene.mad"] <- mad(x)
    res[i,"gene.min"] <- min(x)
    res[i,"gene.max"] <- max(x)
    res[i,"gene.skewness"] <- skewness(x)
    res[i,"gene.kurtosis"] <- kurtosis(x)
    if(length(x)>5) plot(density(x),main=paste("l2fc\n",parent),cex.lab=1.2,cex.axis=1.2,xlab="l2fc",xlim=c(-3,3),cex.main=0.8)
    else plot(c(1,1),main="missing plot")

    temp2 <- abs(temp2)
    cm <- colMeans(temp2)

    genes1 <- names(cm)
    genes2 <- temp1$gene_symbol
    genes <- genes1[is.element(genes1,genes2)]
    temp1 <- temp1[genes,]
    cm <- cm[genes]
    #browser()
    #cat(length(cm),length(temp1$log2_CPM_mean),"\n")
    if(length(cm)!=length(temp1$log2_CPM_mean)) {
      cat("length mismatch\n")
      browser()
    }
    x <- cm/temp1$log2_CPM_mean
    res[i,"l2fccpm.mean"] <- mean(x)
    res[i,"l2fccpm.sd"] <- sd(x)
    res[i,"l2fccpm.med"] <- median(x)
    res[i,"l2fccpm.mad"] <- mad(x)
    res[i,"l2fccpm.min"] <- min(x)
    res[i,"l2fccpm.max"] <- max(x)
    res[i,"l2fccpm.skewness"] <- skewness(x)
    res[i,"l2fccpm.kurtosis"] <- kurtosis(x)
    if(length(x)>5) plot(density(x),main=paste("l2fc/cpm\n",parent),cex.lab=1.2,cex.axis=1.2,xlab="l2fc/cpm",xlim=c(0,0.5),cex.main=0.8)
    else plot(c(1,1),main="missing plot")

    temp <- signaturescoremat[is.element(signaturescoremat$parent,parent),]
    x <- temp$signature_score
    res[i,"score.mean"] <- mean(x)
    res[i,"score.sd"] <- sd(x)
    res[i,"score.med"] <- median(x)
    res[i,"score.mad"] <- mad(x)
    res[i,"score.min"] <- min(x)
    res[i,"score.max"] <- max(x)
    res[i,"score.skewness"] <- skewness(x)
    res[i,"score.kurtosis"] <- kurtosis(x)

    if(length(x)>5) plot(density(x),main=paste("sigscore\n",parent),cex.lab=1.2,cex.axis=1.2,xlab="signature score",xlim=c(-1,1),cex.main=0.8)
    else plot(c(1,1),main="missing plot")
    temp <- signaturecr[is.element(signaturecr$signature,parent),]
    temp <- temp[temp$hitcall>0.5,]
    bmd <- temp$bmd
    res[i,"nhit"] <- length(bmd)
    res[i,"nhit.lt.1"] <- length(bmd[bmd<0.1])
    res[i,"nhit.lt.0.1"] <- length(bmd[bmd<0.01])

    if(!to.file) browser()
    if(i%%100==0) cat("finished ",i,"out of ",np,"\n")
  }
  file <- paste0("../output/signature_corr/signatureScoreVariance.xlsx")
  write.xlsx(res,file)

  if(!to.file) browser()
  else dev.off()
}
