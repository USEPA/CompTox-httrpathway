#--------------------------------------------------------------------------------------
#' Find signatures out of a set of chemicals
#' After running this, run the following ...
#' superTargetPODplot
#' superTargetStats
#'
#' @param to.file If TRUE, send the plots to a file
#' @param do.load If TRUE, load hte large HTTr data set into memory
#' @param dataset Name of the HTTr data set
#' @param celltype Name of the cell type
#' @param ngene the number of genes to consider
#' @param cutoff The minimum number of signatures hat have to be active in a super
#' target for the super target to be considered active. Default is 5
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics par
#' @importFrom openxlsx read.xlsx
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gplots heatmap.2
#' @importFrom stats hclust dist
#'
# dataset="mcf7_ph1_pe1_normal_block_123_allPG",
# sigcatalog="signatureDB_master_catalog 2021-08-27",
# sigset="screen_large",
#
#
# sigcatalog="signatureDB_master_catalog ER",sigset="estrogen"
#' @export signatureFinder
#
#-------------------------------------------------------------------------------
signatureFinder <- function(to.file=F,
                            do.load=F,
                            dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                            celltype="MCF7",
                            ngene=200,
                            cutoff = 0.9,
                            chemfile="../ERModel/ER_chems all mcf7_ph1_pe1_normal_block_123_allPG screen_large 0.9 10.xlsx") {
  printCurrentFunction(paste(dataset))
  if(to.file) {
    fname <- paste0("../input/signatureFinder/signatureFinder ",dataset,".pdf")
    pdf(file=fname,width=8,height=12,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,5,5,2))

  ############################################################################################################
  # read the data
  ############################################################################################################
  if(do.load) {
    cat("Load the data\n")
    file = paste0("../input/fcdata/CHEM_DICT_",dataset,".RDS")
    print(file)
    CHEM_DICT <<- readRDS(file)
    file = paste0("../input/fcdata/FCMAT2_",dataset,".RDS")
    print(file)
    FCMAT2 <<- readRDS(file)
  }

  chems = read.xlsx(chemfile)
  rownames(chems)=chems$dtxsid
  chems[is.na(chems$auc.agonist),"auc.agonist"] = 0
  chems[is.na(chems$auc.antagonist),"auc.antagonist"] = 0
  d1 = chems[chems$auc.agonist>0.4,"dtxsid"]
  d2 = chems[chems$auc.antagonist>0.4,"dtxsid"]

  name.list = c("signature","parent","source","subsource","type","direction","ngene","description",
                "target_class_0","super_target","effect_direction","target_class","super_target_level","sample_key")
  catrow = as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  catrow[1,"source"] = "CCTE"
  catrow[1,"subsource"] = dataset
  catrow[1,"type"] = "nondirectional"
  catrow[1,"direction"] = "nondirectional"
  catrow[1,"description"] = "-"
  catrow[1,"target_class_0"] = "-"
  catrow[1,"super_target"] = "Estrogen"
  catrow[1,"effect_direction"] = "-"
  catrow[1,"target_class"] = "Transcription Factor"
  catrow[1,"super_target_level"] = "gene"
  catrow[1,"ngene"] = ngene

  names(catrow) = name.list

  ####################################################################################
  # agonists
  ####################################################################################
  genelists.up = list()
  genelists.dn = list()
  catalog.up = NULL
  catalog.dn = NULL
  for(dtxsid in d1) {
    logpod = chems[dtxsid,"httr.pod"]
    chem_dict= CHEM_DICT[is.element(CHEM_DICT$dtxsid,dtxsid),]
    logconc = log10(chem_dict[,"conc"])
    useme = -1
    if(logpod<logconc[1]) useme = 1
    else {
      i = which.min(abs(logpod-logconc))
      if(logpod>logconc[i]) useme = i+2
      else useme = i+1
    }

    sk = chem_dict[useme,"sample_key"]
    sig0 = paste(chem_dict[useme,"name"],chem_dict[useme,"dtxsid"],chem_dict[useme,"sample_id"],chem_dict[useme,"conc"],ngene)

    temp = as.data.frame(FCMAT2[sk,])
    temp[,"gene"] = rownames(temp)
    temp = temp[order(temp[,1]),]
    gdn = temp[1:ngene,2]
    sig = paste(sig0,"dn")
    catrow[1,"signature"] = sig
    catrow[1,"parent"] = sig
    catrow[1,"sample_key"] = sk
    catalog.dn = rbind(catalog.dn,catrow)
    genelists.dn[[sig]] = gdn

    temp = temp[order(temp[,1],decreasing=F),]
    gup = temp[1:ngene,2]
    sig = paste(sig0,"up")
    catrow[1,"signature"] = sig
    catrow[1,"parent"] = sig
    catalog.up = rbind(catalog.up,catrow)
    genelists.up[[sig]] = gup
  }

  #######################################################################################
  # plot the l2fc
  #######################################################################################
  sk = catalog.up[,"sample_key"]
  sig = catalog.up[,"signature"]
  fc = t(FCMAT2[sk,])
  colnames(fc) = sig
  fc[is.na(fc)] = 0
  temp = fc
  temp[abs(temp)<0.5] = 0
  temp[abs(temp)>0] = 1
  rs = rowSums(temp)
  fc = fc[rs>1,]
  print(dim(fc))
  fc[fc>1] = 1
  fc[fc< -1] = -1
  heatmap.2(as.matrix(t(fc)),
            margins=c(5,25),
            dendrogram="both",
            scale="none",
            main="l2fc",
            xlab="",
            ylab="",
            cexCol=0.1,
            cexRow=0.8,
            Rowv=T,
            Colv=T,
            trace="none",
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            key=T,
            col=brewer.pal(7,"RdGy"),
            key.title="Key",
            key.xlab="gene",
            cex.main=1)
  if(!to.file) browser()
  #######################################################################################
  # plot the genes
  #######################################################################################
  sigs = catalog.up[,"signature"]
  genes = NULL
  for(sig in sigs) {
    gtemp = genelists.up[sig][[1]]
    genes = c(genes,gtemp)
  }
  gt = table(genes)
  cat("genes:",length(genes),"\n")
  x = gt[gt>=2]
  genes = names(x)
  genes = unique(genes)
  cat("genes:",length(genes),"\n")

  res = as.data.frame(matrix(nrow=length(sigs),ncol=length(genes)))
  rownames(res) = sigs
  names(res) = genes
  res[] = 0
  for(sig in sigs) {
    gtemp = genelists.up[sig][[1]]
    gtemp = gtemp[is.element(gtemp,genes)]
    res[sig,gtemp] = 1
  }
  cs = colSums(res)
  res = res[,cs>1]

  heatmap.2(as.matrix(res),
            margins=c(10,20),
            dendrogram="both",
            scale="none",
            main="agonist up",
            xlab="",
            ylab="",
            cexCol=0.5,
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

  ####################################################################################
  # antagonists
  ####################################################################################
  genelists.up = list()
  genelists.dn = list()
  catalog.up = NULL
  catalog.dn = NULL
  for(dtxsid in d2) {
    logpod = chems[dtxsid,"httr.pod"]
    chem_dict= CHEM_DICT[is.element(CHEM_DICT$dtxsid,dtxsid),]
    logconc = log10(chem_dict[,"conc"])
    useme = -1
    if(logpod<logconc[1]) useme = 1
    else {
      i = which.min(abs(logpod-logconc))
      if(logpod>logconc[i]) useme = i+2
      else useme = i+1
    }

    sk = chem_dict[useme,"sample_key"]
    sig0 = paste(chem_dict[useme,"name"],chem_dict[useme,"dtxsid"],chem_dict[useme,"sample_id"],chem_dict[useme,"conc"],ngene)

    temp = as.data.frame(FCMAT2[sk,])
    temp[,"gene"] = rownames(temp)
    temp = temp[order(temp[,1]),]
    gdn = temp[1:ngene,2]
    sig = paste(sig0,"dn")
    catrow[1,"signature"] = sig
    catrow[1,"parent"] = sig
    catalog.dn = rbind(catalog.dn,catrow)
    genelists.dn[[sig]] = gdn

    temp = temp[order(temp[,1],decreasing=F),]
    gup = temp[1:ngene,2]
    sig = paste(sig0,"up")
    catrow[1,"signature"] = sig
    catrow[1,"parent"] = sig
    catrow[1,"sample_key"] = sk
    catalog.up = rbind(catalog.up,catrow)
    genelists.up[[sig]] = gup
  }

  sigs = catalog.dn[,"signature"]
  genes = NULL
  for(sig in sigs) {
    gtemp = genelists.dn[sig][[1]]
    genes = c(genes,gtemp)
  }
  gt = table(genes)
  cat("genes:",length(genes),"\n")
  x = gt[gt>=1]
  genes = names(x)
  genes = unique(genes)
  cat("genes:",length(genes),"\n")

  res = as.data.frame(matrix(nrow=length(sigs),ncol=length(genes)))
  rownames(res) = sigs
  names(res) = genes
  res[] = 0
  for(sig in sigs) {
    gtemp = genelists.dn[sig][[1]]
    gtemp = gtemp[is.element(gtemp,genes)]
    res[sig,gtemp] = 1
  }
  cs = colSums(res)
  res = res[,cs>1]

  heatmap.2(as.matrix(res),
            margins=c(10,20),
            dendrogram="both",
            scale="none",
            main="antagonist dn",
            xlab="",
            ylab="",
            cexCol=0.5,
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
