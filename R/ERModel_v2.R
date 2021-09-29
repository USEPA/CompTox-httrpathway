#--------------------------------------------------------------------------------------
#' Analyze teh ER model data relative to the MCF7 HTTr data
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
ERModel <- function(to.file=T,
                    do.load=T,
                    dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                    sigcatalog="signatureDB_master_catalog 2021-08-27",
                    sigset="screen_large",
                    method="gsea",
                    celltype="MCF7",
                    hccut=0.9,
                    tccut=1) {
  printCurrentFunction(paste(dataset,sigset,method))

  ############################################################################################################
  # read the data
  ############################################################################################################
  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    mat = mat[is.element(mat$super_target,"Estrogen"),]
    ERSIG <<- mat

    file = paste0("../output/gene_conc_resp_summary/GENE_CR_",dataset,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    GENE <<- GENE_CR

    file = paste0("../input/signatures/signatureDB_genelists.RData")
    print(file)
    load(file=file)
    genelists <<- genelists

    catalog = read.xlsx(paste0("../input/signatures/",sigcatalog,".xlsx"))
    catalog = catalog[catalog[,sigset]==1,]
    catalog <<- catalog
  }
  file = "../ERModel/ER SuperMatrix 2021-09-19 with dtxsid.xlsx"
  supermat = read.xlsx(file)
  supermat = supermat[order(supermat[,"AUC.Agonist"],decreasing=T),]
  supermat$newname = supermat$name

  ############################################################################################################
  # add summary variables
  ############################################################################################################
  ersig = ERSIG
  ersig[is.na(ersig$bmd),"bmd"] = 1000
  ersig[ersig$bmd>1000,"bmd"] = 1000
  ersig$scaled_top = ersig$top_over_cutoff * ersig$top / abs(ersig$top)
  ersig$logbmd = 3-log10(ersig$bmd)
  ersig$auc = ersig$logbmd * ersig$scaled_top
  ersig0 = ersig
  ersig = ersig[ersig$hitcall>=hccut,]
  ersig = ersig[ersig$top_over_cutoff>tccut,]
  erchems = unique(ersig[,c("dtxsid","name")])
  rownames(erchems) = erchems$dtxsid
  for(i in 1:nrow(supermat)) {
    dtxsid = supermat[i,"dtxsid"]
    if(is.element(dtxsid,erchems$dtxsid)) supermat[i,"newname"] = erchems[dtxsid,"name"]
  }

  genecr = GENE
  genecr[is.na(genecr$bmd),"bmd"] = 1000
  genecr[genecr$bmd>1000,"bmd"] = 1000
  genecr$scaled_top = genecr$top_over_cutoff * genecr$top / abs(genecr$top)
  genecr$logbmd = 3-log10(genecr$bmd)
  genecr$auc = genecr$logbmd * genecr$scaled_top
  genecr0 = genecr
  genecr = genecr[genecr$hitcall>=hccut,]
  genecr = genecr[genecr$top_over_cutoff>tccut,]

  dlist1 = unique(supermat$dtxsid)
  dlist2 = unique(ersig$dtxsid)
  cat("ToxCast chemicals:",length(dlist1),"\n")
  cat("HTTr chemicals:",length(dlist2),"\n")

  if(to.file) {
    fname <- paste0("../ERModel/ERmodel.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,15,6,2))
  dlist = supermat[supermat[,"AUC.Agonist"]>=0.4,"dtxsid"]
  dlist = c(dlist,supermat[supermat[,"AUC.Antagonist"]>=0.4,"dtxsid"])
  dlist = dlist[is.element(dlist,ersig$dtxsid)]
  dlist = dlist[!is.element(dlist,"DTXSID8022325")]

  ################################################################################################
  # signature selection
  ################################################################################################
  temp = ersig[is.element(ersig$dtxsid,dlist),]
  chems = unique(temp[,c("dtxsid","name")])
  rownames(chems) = chems$dtxsid
  chems = chems[dlist,]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="logbmd",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat = mat[,chems$name]

  colcol = colnames(mat)
  auclist = vector(mode="integer",length=length(colcol))
  auclist[] = 0
  for(i in 1:length(colcol)) {
    name = colcol[i]
    auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
    auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
    auclist[i] = max(auc1,auc2)
    if(auc2>auc1) auclist[i] = -auclist[i]
  }
  mypal = colorRampPalette( c( "blue", "red" ) )( 10 )
  x = map2color(auclist,mypal)
  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,10),
                      dendrogram="row",
                      scale="none",
                      main="Reference Chemicals\nSignatures BMD",
                      xlab="",
                      ylab="",
                      cexCol=0.8,
                      cexRow=0.1,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="top directieon",
                      cex.main=1,
                      ColSideColors=x)
  #if(!to.file) browser()

  imat = mat
  imat[imat>0] = 1
  rs = rowSums(imat)
  nmin = ncol(imat)
  mask = imat[,1]
  mask[] = 0
  mask[rs>=10] = 1
  sigs = rownames(imat)[mask==1]

  exclude.list = c("CMAP fulvestrant 1e-06 100 5392 100",
                   "CMAP fulvestrant 1e-06 100 6766 100",
                   "CMAP fulvestrant 1e-06 100 8047 100",
                   "CMAP fulvestrant 1e-06 100 8137 100",
                   "CMAP fulvestrant 1e-06 100 8530 100",
                   "CMAP fulvestrant 1e-06 100 9169 100",
                   "CMAP fulvestrant 1e-06 100 9241 100",
                   "CMAP fulvestrant 1e-08 100 1417 100",
                   "CMAP fulvestrant 1e-08 100 5485 100",
                   "CMAP fulvestrant 1e-08 100 8227 100",
                   "CMAP dienestrol 1.5e-05 100 7219 100",
                   "CMAP hexestrol 1.48e-05 100 7342 100",
                   "CMAP diethylstilbestrol 1.5e-05 100 4489 100",
                   "CMAP estradiol 1.46e-05 100 2017 100",
                   "CMAP butyl hydroxybenzoate 2.06e-05 100 6496 100"
                   )
  sigs = sigs[!is.element(sigs,exclude.list)]
  mat = mat[sigs,]

  colcol = colnames(mat)
  auclist = vector(mode="integer",length=length(colcol))
  auclist[] = 0
  for(i in 1:length(colcol)) {
    name = colcol[i]
    auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
    auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
    auclist[i] = max(auc1,auc2)
    if(auc2>auc1) auclist[i]= -auclist[i]
  }
  mypal = colorRampPalette( c( "blue", "red" ) )(10)
  x = map2color(auclist,mypal)

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,20),
                      dendrogram="row",
                      scale="none",
                      main="Reference Chemicals\nBMD Signatures Filtered",
                      xlab="",
                      ylab="",
                      cexCol=0.8,
                      cexRow=0.8,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="log(bmd)",
                      cex.main=1,
                      ColSideColors=x)
  #if(!to.file) browser()

  ################################################################################################
  # scaled_top for the selected signatures and reference chemicals
  ################################################################################################
  temp = ersig0[is.element(ersig0$dtxsid,dlist),]
  temp = temp[is.element(temp$signature,sigs),]
  chems = unique(temp[,c("dtxsid","name")])
  rownames(chems) = chems$dtxsid
  chems = chems[dlist,]

  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="scaled_top",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat = mat[,chems$name]
  mat[is.na(mat)] = 0
  cs = colSums(mat)
  mat = mat[,order(cs)]

  colcol = colnames(mat)
  auclist = vector(mode="integer",length=length(colcol))
  auclist[] = 0
  for(i in 1:length(colcol)) {
    name = colcol[i]
    auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
    auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
    auclist[i] = max(auc1,auc2)
    if(auc2>auc1) auclist[i]= -auclist[i]
  }
  mypal <- colorRampPalette( c( "blue", "red" ) )( 10 )
  x = map2color(auclist,mypal)
  mat[mat>3] = 3
  mat[mat< -3] = -3

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,20),
                      dendrogram="row",
                      scale="none",
                      main="Reference Chemical\nnSignatures Scaled Top",
                      xlab="",
                      ylab="",
                      cexCol=0.8,
                      cexRow=0.8,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top directieon",
                      cex.main=1,
                      ColSideColors=x)
  #if(!to.file) browser()

  ################################################################################################
  # auc for the selected signatures and reference chemicals
  ################################################################################################
  temp = ersig0[is.element(ersig0$dtxsid,dlist),]
  temp = temp[is.element(temp$signature,sigs),]
  chems = unique(temp[,c("dtxsid","name")])
  rownames(chems) = chems$dtxsid
  chems = chems[dlist,]

  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="auc",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat = mat[,chems$name]
  mat[is.na(mat)] = 0
  cs = colSums(mat)
  mat = mat[,order(cs)]
  colcol = colnames(mat)
  auclist = vector(mode="integer",length=length(colcol))
  auclist[] = 0
  for(i in 1:length(colcol)) {
    name = colcol[i]
    auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
    auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
    auclist[i] = max(auc1,auc2)
    if(auc2>auc1) auclist[i]= -auclist[i]
  }
  mypal <- colorRampPalette( c( "blue", "red" ) )( 10 )
  x = map2color(auclist,mypal)
  limit = 10
  mat[mat> limit] = limit
  mat[mat< -limit] = -limit

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,20),
                      dendrogram="row",
                      scale="none",
                      main="Reference Chemical\nnSignatures AUC",
                      xlab="",
                      ylab="",
                      cexCol=0.8,
                      cexRow=0.8,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top directieon",
                      cex.main=1,
                      ColSideColors=x)
  #if(!to.file) browser()

  ################################################################################################
  # Get the ER agonists based on the sig list - bmd
  ################################################################################################
  temp = ersig[is.element(ersig$signature,sigs),]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="logbmd",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  imat = mat
  imat[imat>0] = 1
  cs = colSums(imat)
  mask = imat[1,]
  mask[] = 0
  mask[cs>=6] = 1
  mat = mat[,mask==1]
  chem.names = colnames(mat)
  chems = unique(ersig[,c("dtxsid","name")])
  chems = chems[is.element(chems$name,chem.names),]
  chems$nhit = 0
  chems$mean.logbmd = 0
  chems$auc.agonist = NA
  chems$auc.antagonist = NA
  for(i in 1:nrow(chems)) {
    cname = chems[i,"name"]
    chems[i,"nhit"] = cs[cname][[1]]
    chems[i,"mean.logbmd"] = mean(mat[,cname])
    dtxsid = chems[i,"dtxsid"]
    temp = supermat[is.element(supermat$dtxsid,dtxsid),"AUC.Agonist"]
    chems[i,"auc.agonist"] = temp[1]
    temp = supermat[is.element(supermat$dtxsid,dtxsid),"AUC.Antagonist"]
    chems[i,"auc.antagonist"] = temp[1]
  }
  chems = chems[order(chems$mean.logbmd,decreasing=T),]

  colcol = colnames(mat)
  auclist = vector(mode="integer",length=length(colcol))
  auclist[] = -1
  auclist = as.numeric(auclist)
  for(i in 1:length(colcol)) {
    name = colcol[i]
    if(is.element(name,supermat$newname)) {
      auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
      auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
      auclist[i] = max(auc1,auc2)
      if(auc2>auc1) auclist[i]= -auclist[i]
    }
  }
  mypal = colorRampPalette( c( "blue", "red" ) )( 10 )
  x = map2color(auclist,mypal)
  x[is.element(x,"#0000FF")] = "#FFFFFF"

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,15),
                      dendrogram="both",
                      scale="none",
                      main="ER Candidates\nnSignatures BMD",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.8,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=x)
  file = "../ERModel/ER_chems.xlsx"
  write.xlsx(chems,file)
  #if(!to.file) browser()

  ################################################################################################
  # Get the ER actives based on the sig list - scaled_top
  ################################################################################################
  dlist = chems$dtxsid
  temp = ersig0[is.element(ersig0$signature,sigs),]
  temp = temp[is.element(temp$dtxsid,dlist),]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="scaled_top",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat[is.na(mat)] = 0
  cs = colSums(mat)
  mat = mat[,order(cs)]

  colcol = colnames(mat)
  auclist = vector(mode="integer",length=length(colcol))
  auclist[] = -1
  auclist = as.numeric(auclist)
  for(i in 1:length(colcol)) {
    name = colcol[i]
    if(is.element(name,supermat$newname)) {
      auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
      auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
      auclist[i] = max(auc1,auc2)
      if(auc2>auc1) auclist[i]= -auclist[i]
    }
  }
  mypal = colorRampPalette( c( "blue", "red" ) )( 10 )
  x = map2color(auclist,mypal)
  x[is.element(x,"#0000FF")] = "#FFFFFF"
  limit = 3
  mat[mat>limit] = limit
  mat[mat< -limit] = -limit

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,15),
                      dendrogram="both",
                      scale="none",
                      main="ER Candidates\nnSignatures Scaled Top",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.8,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=x)

  #if(!to.file) browser()

  ################################################################################################
  # Get the ER actives based on the sig list - auc
  ################################################################################################
  dlist = chems$dtxsid
  temp = ersig0[is.element(ersig0$signature,sigs),]
  temp = temp[is.element(temp$dtxsid,dlist),]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="auc",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat[is.na(mat)] = 0
  cs = colSums(mat)
  mat = mat[,order(cs)]
  tmat = as.data.frame(t(mat))
  tmat = cbind(rownames(tmat),tmat,stringsAsFactors=F)
  tmat = cbind(rownames(tmat),tmat,stringsAsFactors=F)
  names(tmat)[1] = "dtxsid"
  names(tmat)[2] = "name"
  tmat = tmat[chems$name,]
  tmat$dtxsid = chems$dtxsid
  file = "../ERModel/ER_chems_auc.xlsx"
  write.xlsx(tmat,file)

  colcol = colnames(mat)
  auclist = vector(mode="integer",length=length(colcol))
  auclist[] = -1
  auclist = as.numeric(auclist)
  for(i in 1:length(colcol)) {
    name = colcol[i]
    if(is.element(name,supermat$newname)) {
      auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
      auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
      auclist[i] = max(auc1,auc2)
      if(auc2>auc1) auclist[i]= -auclist[i]
    }
  }
  mypal = colorRampPalette( c( "blue", "red" ) )( 10 )
  x = map2color(auclist,mypal)
  x[is.element(x,"#0000FF")] = "#FFFFFF"

  limit = 10
  mat[mat>limit] = limit
  mat[mat< -limit] = -limit

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,15),
                      dendrogram="both",
                      scale="none",
                      main="ER Candidates\nnSignatures AUC",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.8,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=x)

  if(!to.file) browser()

  ################################################################################################
  # Plot ER model vs signature stats
  ################################################################################################
  par(mfrow=c(3,2),mar=c(4,4,4,2))
  temp = chems[!is.na(chems$auc.agonist),]
  temp = temp[!is.na(temp$auc.antagonist),]
  x = temp$auc.agonist
  for(i in 1:nrow(temp)) if(temp[i,"auc.antagonist"]>temp[i,"auc.agonist"]) x[i] = temp[i,"auc.antagonist"]
  y = temp$mean.logbmd
  mod = lm(y~x)
  smod = summary(mod)
  intercept = smod$coefficients[1,1]
  slope = smod$coefficients[2,1]
  p = smod$coefficients[2,4]
  plot(y~x,xlab="AUC(agonist)",ylab="Mean(log(bmd))",cex.lab=1.2,cex.axis=1.2,ylim=c(0,4.25),main=paste("p=",format(p,digits=2)))
  for(i in 1:nrow(temp)) {
    ag = temp[i,"auc.agonist"]
    antag = temp[i,"auc.antagonist"]
    col = "red"
    if(antag>ag) col = "blue"
    points(x[i],y[i],pch=21,bg=col)
  }
  lines(c(0,1.5),c(intercept,intercept+slope*1.5))

  x = temp$auc.agonist
  y = temp$nhit
  mod = lm(y~x)
  smod = summary(mod)
  intercept = smod$coefficients[1,1]
  slope = smod$coefficients[2,1]
  p = smod$coefficients[2,4]
  plot(y~x,xlab="AUC(agonist)",ylab="nHit",cex.lab=1.2,cex.axis=1.2,ylim=c(4,12),main=paste("p=",format(p,digits=2)))
  lines(c(0,1.5),c(intercept,intercept+slope*1.5))
  #if(!to.file) browser()

  print("add plot of individual pods, point by signature type")

  ################################################################################################
  # pull out the interesting genes
  ################################################################################################
  genes = NULL
  for(i in 1:length(sigs)) {
    sig = sigs[i]
    temp = catalog[is.element(catalog$parent,sig),"signature"]
    for(x in temp) {
      y = genelists[x][[1]]
      genes = c(genes,y)
    }
  }
  genes = unique(genes)
  genecr = genecr[is.element(genecr$dtxsid,chems$dtxsid),]
  genecr = genecr[is.element(genecr$gene,genes),]
  genecr = genecr[genecr$hitcall>hccut,]
  mat = reshape2::dcast(genecr,gene~name,fill=0,value.var="logbmd",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  imat = mat
  imat[imat>0] = 1
  rs = rowSums(imat)
  mask = imat[,1]
  mask[] = 0
  mask[rs>=10] = 1
  mat = mat[mask==1,]

  colcol = colnames(mat)
  auclist = vector(mode="integer",length=length(colcol))
  auclist[] = -1
  auclist = as.numeric(auclist)
  for(i in 1:length(colcol)) {
    name = colcol[i]
    if(is.element(name,supermat$newname)) {
      auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
      auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
      auclist[i] = max(auc1,auc2)
      if(auc2>auc1) auclist[i]= -auclist[i]
    }
  }
  mypal = colorRampPalette( c( "blue", "red" ) )( 10 )
  x = map2color(auclist,mypal)
  x[is.element(x,"#0000FF")] = "#FFFFFF"

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,10),
                      dendrogram="both",
                      scale="none",
                      main="ER Candidates\nGene BMD",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.2,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=x)
  if(!to.file) browser()

  genecr = genecr0
  genecr = genecr[genecr$hitcall>0.25,]
  mat = reshape2::dcast(genecr,gene~name,fill=0,value.var="scaled_top",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat[is.na(mat)] = 0
  imat = abs(mat)
  imat[imat>0] = 1
  rs = rowSums(imat)
  mask = imat[,1]
  mask[] = 0
  mask[rs>=10] = 1
  mat = mat[mask==1,]
  cs = colSums(mat)
  mat = mat[,order(cs)]

  limit = 5
  mat[mat>limit] = limit
  mat[mat< -limit] = -limit

  colcol = colnames(mat)
  auclist = vector(mode="integer",length=length(colcol))
  auclist[] = -1
  auclist = as.numeric(auclist)
  for(i in 1:length(colcol)) {
    name = colcol[i]
    if(is.element(name,supermat$newname)) {
      auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
      auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
      auclist[i] = max(auc1,auc2)
      if(auc2>auc1) auclist[i]= -auclist[i]
    }
  }
  mypal = colorRampPalette( c( "blue", "red" ) )( 10 )
  x = map2color(auclist,mypal)
  x[is.element(x,"#0000FF")] = "#FFFFFF"
  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,10),
                      dendrogram="both",
                      scale="none",
                      main="ER Candidates\nGene Scaled Top",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.2,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=x)
  if(!to.file) browser()
  else dev.off()
}
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
