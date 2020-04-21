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
noiseHunter.controls.hm <- function(to.file=F,
                                    dir="../input/httr_mcf7_screen/meanncnt0_5-plateteffect_0-shrinkage_normal_DMEM_6_controls/",
                                    ngene=100,
                                    cutoff.l2fcfc=1){
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/noiseHunter/noiseHunter.controls ",cutoff," heatmap.pdf")
    if(ngene>0) fname <- paste0("../output/noiseHunter/noiseHunter.controls. ",cutoff," heatmap ngene ",ngene,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  prefix <- "httr_mcf7_ph1_"
  type.list <- c("QC_BL_TSA","QC_HBRR-vs-UHRR","refChem_GEN","refChem_SIRO","refChem_TSA")
  type.list <- c("refChem_GEN","refChem_SIRO","refChem_TSA")
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
    res <- reshape2::dcast(mat,probe_id~trt_name,value.var="log2FoldChange",fill=0)
    rownames(res) <- res[,1]
    res <- as.matrix(res[,2:ncol(res)])
    if(ngene>0) res <- res[1:ngene,]

    res[res>cutoff.l2fcfc] <- cutoff.l2fcfc
    res[res< -cutoff.l2fcfc] <- -cutoff.l2fcfc
    result <- heatmap.2(res,
                        margins=c(10,10),
                        dendrogram="row",
                        scale="none",
                        main=paste("plate group:",pg),
                        xlab="",
                        ylab="",
                        cexCol=1,
                        cexRow=0.1,
                        Rowv=T,
                        Colv=F,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(9,"PiYG"),
                        key.title="Key",
                        key.xlab="l2fc",
                        cex.main=1)

    if(!to.file) browser()
  }
  if(to.file) dev.off()
}
