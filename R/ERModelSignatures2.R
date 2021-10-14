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
ERModelSignatures2 <- function(to.file=F,
                              dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                              sigcatalog="signatureDB_master_catalog 2021-09-29",
                              sigset="screen_large",
                              method="gsea",
                              hccut=0.9,
                              minhit=10) {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../ERModel/sigs2/ERmodelSignatures2.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(4,5,5,2))

  ############################################################################################################
  # read the data
  ############################################################################################################
  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  mat0 = read.xlsx(file)
  mat0 = mat0[mat0[,sigset]==1,]
  mat0 = mat0[is.element(mat0$super_target,"Estrogen"),]
  #mat = mat0
  mat = mat0[is.element(mat0$source,"CMAP"),]
  file = paste0("../input/signatures/signatureDB_genelists.RData")
  print(file)
  load(file=file)

  ############################################################################################################
  # make the signature x gene matrix
  ############################################################################################################
  sigs = mat[,"signature"]
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
    gtemp = unique(genelists[sig][[1]])
    res[sig,gtemp] = 1
  }
  cs = colSums(res)
  hist(cs,main="signature per gene")
  csdf = as.data.frame(cs)
  csdf$gene = rownames(csdf)
  names(csdf)[1] = "signatures"
  csdf = csdf[,c("gene","signatures")]
  csdf = csdf[order(csdf$signatures,decreasing=T),]
  file = "../ERModel/sigs2/ER gene x sig count.xlsx"
  write.xlsx(csdf,file)
  csdf = csdf[csdf$signatures>=5,]
  res = res[,csdf$gene]
  rs = rowSums(res)
  hist(rs,main="genes per signature")
  rsdf = as.data.frame(rs)
  rsdf$signature = rownames(rsdf)
  names(rsdf)[1] = "genes"
  rsdf = rsdf[,c("signature","genes")]
  rsdf = rsdf[order(rsdf$genes,decreasing=T),]
  res = res[rsdf$signature,]
  file = "../ERModel/sigs2/ER sig x gene count.xlsx"
  write.xlsx(rsdf,file)

  file = "../ERModel/sigs2/signature x gene matrix.RData"
  save(res,file=file)
  heatmap.2(as.matrix(res),
            margins=c(10,5),
            dendrogram="none",
            scale="none",
            main=paste("ER Signatures ",nrow(res),"\nx Genes ",ncol(res)),
            xlab="",
            ylab="",
            cexCol=0.1,
            cexRow=0.1,
            Rowv=F,
            Colv=F,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=F,
            col=brewer.pal(3,"Reds"),
            cex.main=1)

  if(!to.file) browser()
  file = "../ERModel/sigs2/CMAP ER sig classification.xlsx"
  smap = read.xlsx(file)

  ###########################################################################################
  # agonist_up
  ###########################################################################################
  temp = smap[is.element(smap$class,"agonist"),]
  temp = temp[is.element(temp$direction,"up"),]
  sigs = temp$signature
  tres = res[sigs,]
  cs = colSums(tres)
  csdf = as.data.frame(cs)
  csdf$gene = rownames(csdf)
  names(csdf)[1] = "signatures"
  csdf = csdf[,c("gene","signatures")]
  csdf = csdf[order(csdf$signatures,decreasing=T),]
  csdf = csdf[csdf$signatures>=5,]

  tres = tres[,csdf$gene]
  rs = rowSums(tres)
  rsdf = as.data.frame(rs)
  rsdf$signature = rownames(rsdf)
  names(rsdf)[1] = "genes"
  rsdf = rsdf[,c("signature","genes")]
  rsdf = rsdf[order(rsdf$genes,decreasing=T),]
  tres = tres[rsdf$signature,]
  heatmap.2(as.matrix(tres),
            margins=c(5,15),
            dendrogram="none",
            scale="none",
            main=paste("agonist up ",nrow(tres),"\nx Genes ",ncol(tres)),
            xlab="",
            ylab="",
            cexCol=0.1,
            cexRow=0.6,
            Rowv=F,
            Colv=F,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=F,
            col=brewer.pal(3,"Reds"))


  if(!to.file) browser()
  ###########################################################################################
  # agonist_dn
  ###########################################################################################
  temp = smap[is.element(smap$class,"agonist"),]
  temp = temp[is.element(temp$direction,"dn"),]
  sigs = temp$signature
  tres = res[sigs,]
  cs = colSums(tres)
  csdf = as.data.frame(cs)
  csdf$gene = rownames(csdf)
  names(csdf)[1] = "signatures"
  csdf = csdf[,c("gene","signatures")]
  csdf = csdf[order(csdf$signatures,decreasing=T),]
  csdf = csdf[csdf$signatures>=5,]

  tres = tres[,csdf$gene]
  rs = rowSums(tres)
  rsdf = as.data.frame(rs)
  rsdf$signature = rownames(rsdf)
  names(rsdf)[1] = "genes"
  rsdf = rsdf[,c("signature","genes")]
  rsdf = rsdf[order(rsdf$genes,decreasing=T),]
  tres = tres[rsdf$signature,]
  heatmap.2(as.matrix(tres),
            margins=c(5,15),
            dendrogram="none",
            scale="none",
            main=paste("agonist dn",nrow(tres),"\nx Genes ",ncol(tres)),
            xlab="",
            ylab="",
            cexCol=0.1,
            cexRow=0.6,
            Rowv=F,
            Colv=F,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=F,
            col=brewer.pal(3,"Reds"))


  if(!to.file) browser()
  ###########################################################################################
  # antagonist up
  ###########################################################################################
  temp = smap[is.element(smap$class,"antagonist"),]
  temp = temp[is.element(temp$direction,"up"),]
  sigs = temp$signature
  tres = res[sigs,]
  cs = colSums(tres)
  csdf = as.data.frame(cs)
  csdf$gene = rownames(csdf)
  names(csdf)[1] = "signatures"
  csdf = csdf[,c("gene","signatures")]
  csdf = csdf[order(csdf$signatures,decreasing=T),]
  csdf = csdf[csdf$signatures>=5,]

  tres = tres[,csdf$gene]
  rs = rowSums(tres)
  rsdf = as.data.frame(rs)
  rsdf$signature = rownames(rsdf)
  names(rsdf)[1] = "genes"
  rsdf = rsdf[,c("signature","genes")]
  rsdf = rsdf[order(rsdf$genes,decreasing=T),]
  tres = tres[rsdf$signature,]
  heatmap.2(as.matrix(tres),
            margins=c(5,15),
            dendrogram="none",
            scale="none",
            main=paste("antagonist up",nrow(tres),"\nx Genes ",ncol(tres)),
            xlab="",
            ylab="",
            cexCol=0.1,
            cexRow=0.6,
            Rowv=F,
            Colv=F,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=F,
            col=brewer.pal(3,"Reds"))

  ###########################################################################################
  # antagonist dn
  ###########################################################################################
  temp = smap[is.element(smap$class,"antagonist"),]
  temp = temp[is.element(temp$direction,"dn"),]
  sigs = temp$signature
  tres = res[sigs,]
  cs = colSums(tres)
  csdf = as.data.frame(cs)
  csdf$gene = rownames(csdf)
  names(csdf)[1] = "signatures"
  csdf = csdf[,c("gene","signatures")]
  csdf = csdf[order(csdf$signatures,decreasing=T),]
  csdf = csdf[csdf$signatures>=5,]

  tres = tres[,csdf$gene]
  rs = rowSums(tres)
  rsdf = as.data.frame(rs)
  rsdf$signature = rownames(rsdf)
  names(rsdf)[1] = "genes"
  rsdf = rsdf[,c("signature","genes")]
  rsdf = rsdf[order(rsdf$genes,decreasing=T),]
  tres = tres[rsdf$signature,]
  heatmap.2(as.matrix(tres),
            margins=c(5,15),
            dendrogram="none",
            scale="none",
            main=paste("antagonist dn ",nrow(tres),"\nx Genes ",ncol(tres)),
            xlab="",
            ylab="",
            cexCol=0.1,
            cexRow=0.6,
            Rowv=F,
            Colv=F,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=F,
            col=brewer.pal(3,"Reds"))

  if(!to.file) browser()
  else dev.off()
}
