#--------------------------------------------------------------------------------------
#' Visualize the l2fc matrix for a given chemical and specific signatures
#'
#' @param min.ngene Signatures will only be saved if the number of genes is >= this value
#' @param max.ngene Signatures will only be saved if the number of genes is <= this value
#' @return No output.
#' @export
#--------------------------------------------------------------------------------------
chemicalFCvisualize <- function(to.file=F,
                                do.load=F,
                                basedir="../input/fcdata/",
                                dataset="DMEM_6hr_screen_normal_pe_1",
                                name="Acetaminophen",
                                parent.list=c("CMAP nordihydroguaiaretic acid 1e-06 100 5443 100","CMAP flucytosine 3.1e-05 100 7510 100"),
                                conc.list=c()){
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

  cmap <- CHEM_DICT[is.element(CHEM_DICT$name,name),]
  mat <- FCMAT2[cmap$sample_key,]
  gl1up <- genelists[paste0(parent.list[1]," up")][[1]]
  gl1dn <- genelists[paste0(parent.list[1]," dn")][[1]]
  gl2up <- genelists[paste0(parent.list[2]," up")][[1]]
  gl2dn <- genelists[paste0(parent.list[2]," dn")][[1]]

  mask1up <- as.numeric(is.element(colnames(mat),gl1up))
  mask1dn <- as.numeric(is.element(colnames(mat),gl1dn))
  mask2up <- as.numeric(is.element(colnames(mat),gl2up))
  mask2dn <- as.numeric(is.element(colnames(mat),gl2dn))
  mat2 <- rbind(mat,mask1up,mask1dn,mask2up,mask2dn)

  mask <- mask1up+mask1dn+mask2up+mask2dn
  mat2 <- mat2[,mask>0]

  if(to.file) {
    fname <- paste0("../output/signature_corr/chemicalFCvisualize ",name,"_",parent.list[1],"_",parent.list[2],".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  result <- heatmap.2(as.matrix(t(mat2)),
                      margins=c(10,10),
                      dendrogram="row",
                      scale="none",
                      main="",
                      xlab="",
                      ylab="",
                      cexCol=0.7,
                      cexRow=0.3,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"PRGn"),
                      key.title="Key",
                      key.xlab="l2fc",
                      cex.main=1)
   file <- paste0("../output/signature_corr/chemicalFCvisualize ",name,"_",parent.list[1],"_",parent.list[2],".xlsx")
  write.xlsx(t(mat2),file)
  if(!to.file) browser()
  else dev.off()
}
