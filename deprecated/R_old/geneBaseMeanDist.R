#--------------------------------------------------------------------------------------
#' get the base mean distribution for each gene
#'
#' @param to.file If TRUE, plot to a file
#' @param do.read If TRUE, read the input file into memory
#' @param dataset The name of the dataset
#' @importFrom utils flush.console
#' @importFrom openxlsx write.xlsx
#' @importFrom e1071 skewness kurtosis
#' @importFrom stats sd
#'
#' @return No output.
#' @export geneBaseMeanDist
#--------------------------------------------------------------------------------------
geneBaseMeanDist<- function(to.file=F,
                            do.read=F,
                            dataset="DMEM_6hr_screen_normal_pe_1"
                            ){
  printCurrentFunction()

  if(do.read) {
    file <- paste0("../input/fcdata/FCMAT1_",dataset,".RDS")
    print(file)
    FCMAT1 <<- readRDS(file)
    cat("data loaded\n")
  }
  cat("copy FCMAT1 to mat\n")
  flush.console()
  mat <- FCMAT1
  gene.list <- sort(unique(mat$gene))
  #gene.list <- gene.list[1:10]
  name.list <- c("gene","nrow.basemean","nunique.basemean",
                 "mean.basemean","sd.basemean",
                 "kurtosis.basemean","skewness.basemean",
                 "nrow.l2fc","nunique.l2fc",
                 "mean.l2fc","sd.l2fc",
                 "kurtosis.l2fc","skewness.l2fc"
  )
  res <- as.data.frame(matrix(nrow=length(gene.list),ncol=length(name.list)))
  names(res) <- name.list
  res$gene <- gene.list
  rownames(res) <- gene.list
  counter <- 0
  for(gene in gene.list) {
    temp <- mat[is.element(mat$gene,gene),]
    x <- temp$basemean
    x <- x[!is.na(x)]
    res[gene,"mean.basemean"] <- mean(x,na.rm=T)
    res[gene,"sd.basemean"] <- sd(x,na.rm=T)
    res[gene,"kurtosis.basemean"] <- kurtosis (x,na.rm=T)
    res[gene,"skewness.basemean"] <- skewness(x,na.rm=T)
    res[gene,"nrow.basemean"] <- length(x)
    res[gene,"nunique.basemean"] <- length((unique(x)))

    y <- temp$l2fc
    y <- y[!is.na(y)]
    res[gene,"mean.l2fc"] <- mean(y,na.rm=T)
    res[gene,"sd.l2fc"] <- sd(y,na.rm=T)
    res[gene,"kurtosis.l2fc"] <- kurtosis (y,na.rm=T)
    res[gene,"skewness.l2fc"] <- skewness(y,na.rm=T)
    res[gene,"nrow.l2fc"] <- length(y)
    res[gene,"nunique.l2fc"] <- length((unique(y)))
    counter <- counter+1
    cat(gene,counter,length(gene.list),"\n")
  }
  file <- paste0("../output/signature_corr/geneBaseMeanDist.xlsx")
  write.xlsx(res,file)
}
