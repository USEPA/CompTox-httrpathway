#--------------------------------------------------------------------------------------
#' Heatmap of the low conc rows
#'
#' @param min.ngene Signatures will only be saved if the number of genes is >= this value
#' @param max.ngene Signatures will only be saved if the number of genes is <= this value
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
lowConcHM <- function(to.file=F,
                      do.load=F,
                      basedir="../input/fcdata/",
                      dataset="DMEM_6hr_screen_normal_pe_1",
                      cutoff=1.5,
                      do.heatmap=F){
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
  }
  FCMAT2[is.nan(FCMAT2)] <- 0
  mat <- FCMAT2
  cmap <- CHEM_DICT[CHEM_DICT$conc_index==1,]
  mat <- mat[cmap$sample_key,]
  hit <- mat
  hit[abs(hit)<cutoff] <- 0
  hit[abs(hit)>0] <- 1
  print(dim(hit))
  rs <- abs(rowSums(hit))
  hit <- hit[rs>0,]
  print(dim(hit))
  cs <- abs(colSums(hit))
  hit <- hit[,cs>0]
  print(dim(hit))
  mat <- mat[rownames(hit),colnames(hit)]
  if(to.file) {
    fname <- paste0("../output/signature_corr/lowConcHM.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  if(do.heatmap) {
    mat2 <- mat#[,1:100]
    result <- heatmap.2(as.matrix(t(mat2)),
                        margins=c(10,10),
                        dendrogram="both",
                        scale="none",
                        main="",
                        xlab="",
                        ylab="",
                        cexCol=0.1,
                        cexRow=0.2,
                        Rowv=T,
                        Colv=T,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(9,"PRGn"),
                        key.title="Key",
                        key.xlab="l2fc",
                        cex.main=1)

    if(!to.file) browser()
  }
  par(mfrow=c(5,2),mar=c(2,4,4,4))

  file <- "../input/chemicals/HTTr.MCF7.Screen.Sample.Key.20180605.xlsx"
  sk <- read.xlsx(file)
  for(key in rownames(mat)) {
    chem <- CHEM_DICT[key,"name"]
    sample_id <- CHEM_DICT[key,"sample_id"]
    block <- unique(sk[is.element(sk$EPA_Sample_ID,sample_id),"Block_Number"])
    main <- paste(chem,":",block)
    plot(density(mat[key,]),main=main,xlim=c(-3,3),xlab="")
    if(!to.file) browser()
  }



  #par(mfrow=c(5,4),mar=c(2,1,2,2))
  #for(gene in colnames(mat)) {
  #  plot(density(mat[,gene]),main=gene,xlim=c(-3,3),xlab="")
  #  browser()
  #}
  if(!to.file) browser()
  else dev.off()
}
