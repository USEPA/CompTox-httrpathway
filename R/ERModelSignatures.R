#' Look at correlation of genes between ER signatures
#'
#' @param to.file If TRUE, send the plots to a file
#' @param do.load If TRUE, load hte large HTTr data set into memory
#' @param dataset Name of the HTTr data set
#' @param sigcatalog Name of the signature catalog to use
#' @param sigset Name of the signature set
#' @param method Scoring method
#' @param celltype Name of the cell type
#' @param hccut Exclude rows in the data set with hitcall less than this value
#' @param tccut Exclude rows in the data set with top_over_cutoff less than this value
#' @param cutoff The minimum number of signatures hat have to be active in a super
#' target for the super target to be considered active. Default is 5
#' @param minconc Minimum concentration used in the plots
#' @param maxconc Maximum concentration used in the plots
#'
#' After running this, run the following ...
#' superTargetPODplot
#' superTargetStats
#-------------------------------------------------------------------------------
ERModelSignatures <- function(to.file=F,
                              dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                              sigcatalog="signatureDB_master_catalog ER",
                              sigset="estrogen",
                              method="gsea",
                              hccut=0.9,
                              minhit=10) {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../ERModel/ERmodelSignatures.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,1),mar=c(4,15,6,2))

  file = paste0("../ERModel/ER_gene_diff ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  genediff = read.xlsx(file)
  genediff = genediff[order(genediff$diff,decreasing=T),]
  diffgenes = genediff[genediff$diff>=1,"gene"]
  ############################################################################################################
  # read the data
  ############################################################################################################
  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  mat0 = read.xlsx(file)
  mat = mat0[is.element(mat0$source,"CMAP"),]
  file = paste0("../input/signatures/signatureDB_genelists.RData")
  print(file)
  load(file=file)

  ############################################################################################################
  # agonist_up
  ############################################################################################################
  cat("agonist up\n")
  sigs = mat[mat[,"agonist_up"]==1,"signature"]
  genes = NULL
  for(sig in sigs) {
    gtemp = genelists[sig][[1]]
    genes = c(genes,gtemp)
  }
  genes = unique(genes)
  cat("genes:",length(genes),"\n")
  res = as.data.frame(matrix(nrow=length(sigs),ncol=length(genes)))
  rownames(res) = sigs
  names(res) = genes
  res[] = 0
  for(sig in sigs) {
    gtemp = genelists[sig][[1]]
    res[sig,gtemp] = 1
  }
  cs = colSums(res)
  res = res[,cs>1]
  cs = colSums(res)
  hist(cs,main="agonist up")
  heatmap.2(as.matrix(res),
            margins=c(10,10),
            dendrogram="both",
            scale="none",
            main="agonist up",
            xlab="",
            ylab="",
            cexCol=0.1,
            cexRow=0.5,
            Rowv=T,
            Colv=T,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=T,
            col=brewer.pal(3,"Reds"),
            key.title="Key",
            key.xlab="gene",
            cex.main=1)

  dg = diffgenes[is.element(diffgenes,names(res))]
  res2 = res[,dg]
  heatmap.2(as.matrix(res2),
            margins=c(10,10),
            dendrogram="both",
            scale="none",
            main="agonist up, diff genes",
            xlab="",
            ylab="",
            cexCol=0.8,
            cexRow=0.5,
            Rowv=T,
            Colv=T,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=T,
            col=brewer.pal(3,"Reds"),
            key.title="Key",
            key.xlab="gene",
            cex.main=1)
  if(!to.file) browser()
  ############################################################################################################
  # agonist_dn
  ############################################################################################################
  cat("agonist dn\n")
  sigs = mat[mat[,"agonist_dn"]==1,"signature"]
  genes = NULL
  for(sig in sigs) {
    gtemp = genelists[sig][[1]]
    genes = c(genes,gtemp)
  }
  genes = unique(genes)
  cat("genes:",length(genes),"\n")
  res = as.data.frame(matrix(nrow=length(sigs),ncol=length(genes)))
  rownames(res) = sigs
  names(res) = genes
  res[] = 0
  for(sig in sigs) {
    gtemp = genelists[sig][[1]]
    res[sig,gtemp] = 1
  }
  cs = colSums(res)
  res = res[,cs>1]
  cs = colSums(res)
  hist(cs,main="agonist dn")
  heatmap.2(as.matrix(res),
            margins=c(10,10),
            dendrogram="both",
            scale="none",
            main="agonist dn",
            xlab="",
            ylab="",
            cexCol=0.1,
            cexRow=0.5,
            Rowv=T,
            Colv=T,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=T,
            col=brewer.pal(3,"Reds"),
            key.title="Key",
            key.xlab="gene",
            cex.main=1)
  dg = diffgenes[is.element(diffgenes,names(res))]
  res2 = res[,dg]
  heatmap.2(as.matrix(res2),
            margins=c(10,10),
            dendrogram="both",
            scale="none",
            main="agonist dn, diff genes",
            xlab="",
            ylab="",
            cexCol=0.8,
            cexRow=0.5,
            Rowv=T,
            Colv=T,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=T,
            col=brewer.pal(3,"Reds"),
            key.title="Key",
            key.xlab="gene",
            cex.main=1)
  if(!to.file) browser()
  ############################################################################################################
  # antagonist_up
  ############################################################################################################
  cat("antagonist up\n")
  sigs = mat[mat[,"antagonist_up"]==1,"signature"]
  genes = NULL
  for(sig in sigs) {
    gtemp = genelists[sig][[1]]
    genes = c(genes,gtemp)
  }
  genes = unique(genes)
  cat("genes:",length(genes),"\n")
  res = as.data.frame(matrix(nrow=length(sigs),ncol=length(genes)))
  rownames(res) = sigs
  names(res) = genes
  res[] = 0
  for(sig in sigs) {
    gtemp = genelists[sig][[1]]
    res[sig,gtemp] = 1
  }
  cs = colSums(res)
  res = res[,cs>1]
  cs = colSums(res)
  hist(cs,main="antagonist up")
  heatmap.2(as.matrix(res),
            margins=c(10,10),
            dendrogram="both",
            scale="none",
            main="antagonist up",
            xlab="",
            ylab="",
            cexCol=0.1,
            cexRow=0.5,
            Rowv=T,
            Colv=T,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=T,
            col=brewer.pal(3,"Reds"),
            key.title="Key",
            key.xlab="gene",
            cex.main=1)
  dg = diffgenes[is.element(diffgenes,names(res))]
  res2 = res[,dg]
  heatmap.2(as.matrix(res2),
            margins=c(10,10),
            dendrogram="both",
            scale="none",
            main="antagonist up, diff genes",
            xlab="",
            ylab="",
            cexCol=0.8,
            cexRow=0.5,
            Rowv=T,
            Colv=T,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=T,
            col=brewer.pal(3,"Reds"),
            key.title="Key",
            key.xlab="gene",
            cex.main=1)
  if(!to.file) browser()
  ############################################################################################################
  # antagonist_dn
  ############################################################################################################
  cat("antagonist dn\n")
  sigs = mat[mat[,"antagonist_dn"]==1,"signature"]
  genes = NULL
  for(sig in sigs) {
    gtemp = genelists[sig][[1]]
    genes = c(genes,gtemp)
  }
  genes = unique(genes)
  cat("genes:",length(genes),"\n")
  res = as.data.frame(matrix(nrow=length(sigs),ncol=length(genes)))
  rownames(res) = sigs
  names(res) = genes
  res[] = 0
  for(sig in sigs) {
    gtemp = genelists[sig][[1]]
    res[sig,gtemp] = 1
  }
  cs = colSums(res)
  res = res[,cs>1]
  cs = colSums(res)
  hist(cs,main="antagonist dn")
  heatmap.2(as.matrix(res),
            margins=c(10,10),
            dendrogram="both",
            scale="none",
            main="antagonist dn",
            xlab="",
            ylab="",
            cexCol=0.1,
            cexRow=0.5,
            Rowv=T,
            Colv=T,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=T,
            col=brewer.pal(3,"Reds"),
            key.title="Key",
            key.xlab="gene",
            cex.main=1)
  dg = diffgenes[is.element(diffgenes,names(res))]
  res2 = res[,dg]
  heatmap.2(as.matrix(res2),
            margins=c(10,10),
            dendrogram="both",
            scale="none",
            main="antagonist dn, diff genes",
            xlab="",
            ylab="",
            cexCol=0.8,
            cexRow=0.5,
            Rowv=T,
            Colv=T,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=T,
            col=brewer.pal(3,"Reds"),
            key.title="Key",
            key.xlab="gene",
            cex.main=1)
  if(!to.file) browser()
  else dev.off()
}
