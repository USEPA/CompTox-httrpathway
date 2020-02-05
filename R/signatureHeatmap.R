library(RColorBrewer)
library(gplots)
#--------------------------------------------------------------------------------------
#
#' Create heatmaps of the htpp data frames
#' @param to.file If TRUE, write plots to a file
#' @param dataset The set of data to be included

#-------------------------------------------------------------------------------------
signatureHeatmap <- function(to.file=F,
                           pathset="PathwaySet_20191025",
                           dataset="DMEM_6hr_pilot_normal_00",
                           method="fc",
                           conthits=T,
                           nametag=NULL,
                           pval=0.05) {
  printCurrentFunction()
  if(is.null(nametag) && conthits) nametag = "conthits"
  if(to.file) {
    file <- paste0("../output/signature_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_", pval,"_",nametag,".pdf")
    pdf(file=file,width=7,height=7,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }


  file <- paste0("../output/signature_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_", pval,"_", nametag ,".RData")
  print(file)
  load(file)
  mat <- PATHWAY_CR

  temp <- mat[,c("name","signature","bmdl")]
  x <- temp$bmdl
  x[is.na(x)] <- 1000
  x <- 3-log10(x)
  x[mat$hitcall<0.5] <- 0
  temp$bmdl <- x
  res <- reshape2::dcast(temp,name~signature)
  rownames(res) <- res[,1]
  res <- res[,2:ncol(res)]

  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  cindex <- read.xlsx(file)
  nlist <- rownames(res)
  rownames(cindex) <- cindex$name
  ccolors <- cindex[nlist,"color"]

  file = paste0("../input/processed_signature_data/PATHWAY_CATALOG_",pathset,".RData")
  load(file=file)
  pindex <- signature_catalog
  colors <- names(res)
  colors[] <- "white"
  for(i in 1:length(colors)) {
    signature <- names(res)[i]
    super_class <- pindex[is.element(pindex$signature,signature),"super_class"]
    if(super_class=="estrogen") colors[i] <- "blue"
    if(is.element(super_class,c("apoptosis","nfkb","oxidative stress","mitochondria","cytotoxicity",
                           "hif1a","hypoxia","extracellular matrix","heat shock","microtubule","h1f2a"))) colors[i] <- "red"
    if(super_class=="cell cycle") colors[i] <- "cyan"
    if(super_class=="p450") colors[i] <- "green"
    if(super_class=="steroid synthesis") colors[i] <- "brown"

  }

  res[res>6] <- 6
  result <- heatmap.2(as.matrix(res),
                      margins=c(10,10),
                      scale="none",
                      main=paste0(dataset),
                      cex.main=0.9,
                      xlab="",
                      ylab="",
                      cexCol=0.1,
                      cexRow=0.5,
                      col=brewer.pal(8,"Reds"),
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      key.title="Key",
                      key.xlab="log(AC50)",
                      ColSideColors=colors,
                      RowSideColors=ccolors)


  if(!to.file) browser()
  if(to.file) dev.off()
}
