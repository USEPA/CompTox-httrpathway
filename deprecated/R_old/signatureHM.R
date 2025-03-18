#--------------------------------------------------------------------------------------
#' Visualize the variance by signature - this supplants signatureVariance.R
#'
#' @param min.ngene Signatures will only be saved if the number of genes is >= this value
#' @param max.ngene Signatures will only be saved if the number of genes is <= this value
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
signatureHM <- function(to.file=F,
                        do.load=F,
                        sigset="pilot_large_all_CMAP",
                        sigcatalog="signatureDB_master_catalog 2020-01-31",
                        dataset="DMEM_6hr_screen_normal_pe_1",
                        method="mygsea"){
  printCurrentFunction()

  if(do.load) {
    file <- paste("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData",sep="")
    print(file)
    load(file=file)
    signaturecr <<- SIGNATURE_CR
    annotations <<- signatureCatalogLoader(sigset,sigcatalog)

    file <- "../input/chemicals/HTTr.MCF7.Screen.Sample.Key.20180605.xlsx"
    bmap <- read.xlsx(file)
    bmap <- bmap[,c("EPA_Sample_ID","Block_Number")]
    names(bmap) <- c("sample_id","block")
    bmap <- bmap[!is.na(bmap$sample_id),]
    bmap <- unique(bmap)
    rownames(bmap) <- bmap$sample_id
    bmap <<- bmap

    file <- "../input/chemicals/HTTr_screen_sample_map.xlsx"
    smap <- read.xlsx(file)
    smap <- smap[,c("sample_id","dtxsid","casrn","name")]
    smap <- unique(smap)
    smap <<- smap
  }

  if(to.file) {
    fname <- paste0("../output/signature_corr/signatureHM.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  temp <- annotations[,c("parent","target_class","super_target")]
  temp <- unique(temp)
  rownames(temp) <-temp$parent

  mat <- signaturecr
  #mat <- mat[1:10000,]
  x <- mat$hitcall
  x[x<0.5] <- 0
  x[x>0] <- 1
  mat$hitcall <- x
  res <- reshape2::dcast(mat,name~signature,value.var="hitcall",fun.aggregate = mean,fill=0)
  rownames(res) <- res[,1]
  res <- res[,2:ncol(res)]

  cs <- colSums(res)
  ns <- length(cs)
  res2 <- as.data.frame(matrix(ncol=2,nrow=ns))
  names(res2) <- c("signature","hits")
  res2$signature <- names(cs)
  res2$hits <- cs
  rownames(res2) <- res2$signature
  temp <- temp[res2$signature,]
  res2 <- cbind(res2,temp)

  file <- "../output/signature_corr/signatureHitCounts.xlsx"
  write.xlsx(res2,file)
  cat("file written, start heatmap\n")
  mat <- signaturecr[signaturecr$hitcall>0.5,]
  #mat <- mat[1:10000,]
  mat$block <- bmap[mat$sample_id,"block"]

  res <- reshape2::dcast(mat,sample_id~signature,value.var="bmd",fun.aggregate = mean,fill=1000)
  rownames(res) <- res[,1]
  res <- res[,2:ncol(res)]
  res <- 3-log10(res)
  block <- bmap[rownames(res),"block"]
  colormap <- vector(length=length(block),mode="character")
  colormap[] <- "white"
  colormap[block==1] <- "blue"
  colormap[block==2] <- "cyan"
  colormap[block==3] <- "green"
  colormap[block==5] <- "red"
  #mat <- mat[,1:20]
  result <- heatmap.2(as.matrix(res),
                      margins=c(10,10),
                      scale="none",
                      main=paste(dataset,"\n",method,":",sigset),
                      xlab="",
                      ylab="",
                      cexCol=0.1,
                      cexRow=0.1,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="l2fc",
                      cex.main=1,
                      RowSideColors=colormap)

  if(!to.file) browser()
  else dev.off()
}
