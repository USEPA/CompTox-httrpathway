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
#'
#'
#' After running this, run the following ...
#' superTargetPODplot
#' superTargetStats
# dataset="mcf7_ph1_pe1_normal_block_123_allPG",
# sigcatalog="signatureDB_master_catalog 2021-08-27",
# sigset="screen_large",
#
#
# sigcatalog="signatureDB_master_catalog ER",sigset="estrogen"
#
#
#-------------------------------------------------------------------------------
ERModel2 <- function(to.file=F,
                    do.stop=T,
                    do.load.sigs=F,
                    do.load.genes=F,
                    dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                    sigcatalog="signatureDB_master_catalog 2021-08-27",
                    sigset="screen_large",
                    method="gsea",
                    celltype="MCF7",
                    hccut=0.9,
                    tccut=1,
                    minhit=10) {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../ERModel/ERmodel ",dataset," ",sigset," ",hccut," ",minhit,".pdf")
    pdf(file=fname,width=8,height=8,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,15,6,2))

  ############################################################################################################
  # read the data
  ############################################################################################################
  if(do.load.sigs) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    mat = mat[is.element(mat$super_target,"Estrogen"),]
    ERSIG <<- mat

    file = paste0("../input/signatures/signatureDB_genelists.RData")
    print(file)
    load(file=file)
    genelists <<- genelists

    catalog = read.xlsx(paste0("../input/signatures/",sigcatalog,".xlsx"))
    catalog = catalog[catalog[,sigset]==1,]
    catalog <<- catalog
  }

  if(do.load.genes) {
    file = paste0("../output/gene_conc_resp_summary/GENE_CR_",dataset,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    GENE <<- GENE_CR
  }

  file = "../ERModel/ER SuperMatrix 2021-09-19 with dtxsid.xlsx"
  supermat = read.xlsx(file)
  supermat = supermat[is.element(supermat$dtxsid,ERSIG$dtxsid),]
  supermat = supermat[order(supermat[,"AUC.Antagonist"],decreasing=T),]
  supermat = supermat[order(supermat[,"AUC.Agonist"],decreasing=T),]
  supermat$newname = supermat$name

  file = "../ERModel/ER refchems.xlsx"
  refchems = read.xlsx(file)
  rownames(refchems) = refchems$dtxsid
  refchems = refchems[is.element(refchems$dtxsid,ERSIG$dtxsid),]

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
  erchems = unique(ersig0[,c("dtxsid","name")])
  rownames(erchems) = erchems$dtxsid
  supermat$pod.agonist =  3
  supermat$pod.antagonist = 3
  supermat$pod.agonist.10 =  3
  supermat$pod.antagonist.10 = 3
  supermat$nactive.agonist = 0
  supermat$nactive.antagonist = 0
  supermat$mode = "inactive"
  for(i in 1:nrow(supermat)) {
    dtxsid = supermat[i,"dtxsid"]
    if(is.element(dtxsid,erchems$dtxsid)) supermat[i,"newname"] = erchems[dtxsid,"name"]
    auc1 = supermat[i,"AUC.Agonist"]
    auc2 = supermat[i,"AUC.Antagonist"]
    mode = "agonist"
    if(auc2>auc1) mode= "antagonist"
    if(auc1>=0.1 || auc2>=0.1) supermat[i,"mode"] = mode
    a1 = c(
      "NVS_NR_bER_AC50",
      "NVS_NR_hER_AC50",
      "NVS_NR_mERa_AC50",
      "OT_ER_ERaERa_0480_AC50",
      "OT_ER_ERaERa_1440_AC50",
      "OT_ER_ERaERb_0480_AC50",
      "OT_ER_ERaERb_1440_AC50",
      "OT_ER_ERbERb_0480_AC50",
      "OT_ER_ERbERb_1440_AC50",
      "OT_ERa_EREGFP_0120_AC50",
      "OT_ERa_EREGFP_0480_AC50",
      "ATG_ERa_TRANS_up_AC50",
      "ATG_ERE_CIS_up_AC50",
      "Tox21_ERa_BLA_Agonist_ratio_AC50",
      "Tox21_ERa_LUC_BG1_Agonist_AC50",
      "ACEA_T47D_80hr_Positive_AC50"
    )
    a2 = c(
      "NVS_NR_bER_AC50",
      "NVS_NR_hER_AC50",
      "NVS_NR_mERa_AC50",
      "OT_ER_ERaERa_0480_AC50",
      "OT_ER_ERaERa_1440_AC50",
      "OT_ER_ERaERb_0480_AC50",
      "OT_ER_ERaERb_1440_AC50",
      "OT_ER_ERbERb_0480_AC50",
      "OT_ER_ERbERb_1440_AC50",
      "OT_ERa_EREGFP_0120_AC50",
      "OT_ERa_EREGFP_0480_AC50",
      "Tox21_ERa_BLA_Agonist_ratio_AC50",
      "Tox21_ERa_LUC_BG1_Agonist_AC50",
      "Tox21_ERa_BLA_Antagonist_ratio_AC50",
      "Tox21_ERa_LUC_BG1_Antagonist_AC50"
    )
    if(mode=="agonist") {
      temp = supermat[i,a1]
      temp = temp[temp<1000000]
      supermat[i,"nactive.agonist"] = length(temp)
      if(length(temp)>=2) {
        temp = log10(temp)
        supermat[i,"pod.agonist"] = median(temp)
        q = quantile(temp,probs=seq(0,1,0.1))
        supermat[i,"pod.agonist.10"] = q[2]
      }
    }
    if(mode=="antagonist") {
      temp = supermat[i,a2]
      temp = temp[temp<1000000]
      supermat[i,"nactive.antagonist"] = length(temp)
      if(length(temp)>=2) {
        temp = log10(temp)
        supermat[i,"pod.antagonist"] = median(temp)
        q = quantile(temp,probs=seq(0,1,0.1))
        supermat[i,"pod.antagonist.10"] = q[2]
      }
    }
  }
  par(mfrow=c(1,1),mar=c(4,4,4,2))
  x = supermat$AUC.Agonist
  y = supermat$pod.agonist
  z = supermat$nactive.agonist
  x = x[z>=5]
  y = y[z>=5]
  plot(y~x,xlab="AUC",ylab="Mean(log(bmd uM))",cex.lab=1.2,cex.axis=1.2,ylim=c(-3,3),main="ER Model Results",pch=21,bg="blue")
  x = supermat$AUC.Antagonist
  y = supermat$pod.antagonist
  z = supermat$nactive.antagonist
  x = x[z>=5]
  y = y[z>=5]
  points(y~x,pch=21,bg="red")
  if(!to.file && do.stop) browser()

  #######################################################################################################
  # prep GENECR
  #######################################################################################################
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

  dlist = supermat[supermat[,"AUC.Agonist"]>=0.4,"dtxsid"]
  dlist = c(dlist,supermat[supermat[,"AUC.Antagonist"]>=0.4,"dtxsid"])
  dlist = dlist[is.element(dlist,ersig$dtxsid)]
  dlist = dlist[!is.element(dlist,"DTXSID8022325")]

  #######################################################################################################
  # build the erchem matrix - some filed not yet filled in
  #######################################################################################################
  erchems = unique(ersig0[,c("dtxsid","name")])
  erchems$nhit = 0
  erchems$nhit.agonist = 0
  erchems$nhit.antagonist = 0
  erchems$httr.pod = 0
  erchems$httr.mode = NA
  erchems$auc.agonist = NA
  erchems$auc.antagonist = NA
  erchems$hts.mode = NA
  erchems$hts.pod.agonist = NA
  erchems$hts.pod.antagonist = NA
  erchems$hts.pod.agonist.10 = NA
  erchems$hts.pod.antagonist.10 = NA
  erchems$hts.nactive.agonist = NA
  erchems$hts.nactive.antagonist = NA
  erchems$refchem.invitro.agonist = NA
  erchems$refchem.invivo.agonist = NA
  erchems$refchem.invitro.antagonist = NA
  erchems$candidate = 0
  for(i in 1:nrow(erchems)) {
    dtxsid = erchems[i,"dtxsid"]
    if(is.element(dtxsid,supermat$dtxsid)) {
      temp = supermat[is.element(supermat$dtxsid,dtxsid),"AUC.Agonist"]
      erchems[i,"auc.agonist"] = max(temp)
      temp = supermat[is.element(supermat$dtxsid,dtxsid),"AUC.Antagonist"]
      erchems[i,"auc.antagonist"] = max(temp)

      temp = supermat[is.element(supermat$dtxsid,dtxsid),"pod.agonist"]
      erchems[i,"hts.pod.agonist"] = min(temp)
      temp = supermat[is.element(supermat$dtxsid,dtxsid),"pod.antagonist"]
      erchems[i,"hts.pod.antagonist"] = min(temp)

      temp = supermat[is.element(supermat$dtxsid,dtxsid),"pod.agonist.10"]
      erchems[i,"hts.pod.agonist.10"] = min(temp)
      temp = supermat[is.element(supermat$dtxsid,dtxsid),"pod.antagonist.10"]
      erchems[i,"hts.pod.antagonist.10"] = min(temp)

      temp = supermat[is.element(supermat$dtxsid,dtxsid),"nactive.agonist"]
      erchems[i,"hts.nactive.agonist"] = max(temp)
      temp = supermat[is.element(supermat$dtxsid,dtxsid),"nactive.antagonist"]
      erchems[i,"hts.nactive.antagonist"] = max(temp)

      temp = supermat[is.element(supermat$dtxsid,dtxsid),"mode"]
      erchems[i,"hts.mode"] = temp[1]
    }

    if(is.element(dtxsid,refchems$dtxsid)) {
      erchems[i,"refchem.invitro.agonist"] = refchems[dtxsid,"invitro.agonist"]
      erchems[i,"refchem.invivo.agonist"] = refchems[dtxsid,"invivo.agonist"]
      erchems[i,"refchem.invitro.antagonist"] = refchems[dtxsid,"invitro.antagonist"]
    }
  }

  x = erchems$auc.agonist
  y = erchems$auc.antagonist
  x[is.na(x)] = 0
  y[is.na(y)] = 0
  z = x-y
  erchems$order = z
  erchems = erchems[order(erchems$order,decreasing=T),]

  ################################################################################################
  # signature selection - version A
  ################################################################################################
  temp = ersig0[is.element(ersig0$dtxsid,dlist),]
  temp = temp[temp$bmd<10,]
  chems = unique(temp[,c("dtxsid","name")])
  rownames(chems) = chems$dtxsid
  chems = chems[dlist,]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="hitcall",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat = mat[,chems$name]
  mat[is.na(mat)] = 0
  mat[mat<hccut] = 0
  mat[mat>0] = 1
  names = colnames(mat)
  col.colors = names
  col.colors[] = "white"
  for(i in 1:length(names)) {
    name = names[i]
    if(is.element(name,supermat$newname)) {
      auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
      auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
      if(auc1>0.1 && auc2<0.1) col.colors[i] = "blue"
      else if(auc2>0.1 && auc1<0.1) col.colors[i] = "red"
      else if(auc2>0.1 && auc1>0.1) {
        if(auc1>auc2) col.colors[i] = "blue"
        else col.colors[i] = "red"
      }
      else col.colors[i] = "gray"
    }
  }

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,20),
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
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()

  ################################################################################################
  # signature selection
  ################################################################################################

  rs = rowSums(mat)
  hist(rs,main="Chemicals (out of 25) active for each signature")
  # nmin = ncol(imat)
  mask = mat[,1]
  mask[] = 0
  mask[rs>=minhit] = 1
  sigs = rownames(mat)[mask==1]

  exclude.list = c("CMAP butyl hydroxybenzoate 2.06e-05 100 5518 100",
                   "CMAP butyl hydroxybenzoate 2.06e-05 100 6496 100",
                   "CMAP genistein 1e-05 100 475 100",
                   "CMAP fulvestrant 1e-06 100 5392 100",
                   "CMAP fulvestrant 1e-06 100 6766 100",
                   "CMAP fulvestrant 1e-06 100 8047 100",
                   "CMAP fulvestrant 1e-06 100 8137 100",
                   "CMAP fulvestrant 1e-06 100 8530 100",
                   "CMAP fulvestrant 1e-06 100 9169 100",
                   "CMAP fulvestrant 1e-06 100 9241 100",
                   "CMAP fulvestrant 1e-08 100 1417 100",
                   "CMAP fulvestrant 1e-08 100 5485 100",
                   "CMAP fulvestrant 1e-08 100 8227 100"
  )
  exclude.list = c("")
  sigs = sigs[!is.element(sigs,exclude.list)]
  if(sigset=="screen_large" && sigcatalog=="signatureDB_master_catalog 2021-08-27") {
    index = grep("fulvestrant",sigs)
  }
  else if(sigset=="estrogen" && sigcatalog=="signatureDB_master_catalog ER") {
    index1 = grep("fulvestrant",sigs)
    index2 = grep("tamoxifen",sigs)
    index3 = grep("raloxifene",sigs)
    index4 = grep("epitiostanol",sigs)
    index = c(index1,index2,index3,index3)
  }
  sigs.antagonist = sigs[index]
  sigs.agonist = sigs[!is.element(sigs,sigs.antagonist)]

  #######################################################################################################
  # Add the next pieces to erchems - nhit, httr.pod, httr.mode
  #######################################################################################################
  temp = ersig[is.element(ersig$signature,sigs),]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="logbmd",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat1 = mat[sigs.agonist,]
  imat1 = abs(mat1)
  imat1[imat1>0] = 1
  cs1 = colSums(imat1)

  mat2 = mat[sigs.antagonist,]
  imat2 = abs(mat2)
  imat2[imat2>0] = 1
  cs2 = colSums(imat2)
  cat("selected chems 2:",ncol(mat),"\n")

  for(i in 1:nrow(erchems)) {
    dtxsid = erchems[i,"dtxsid"]
    cname = erchems[i,"name"]
    erchems[i,"httr.pod"] = 3
    if(is.element(cname,names(mat))) {
      erchems[i,"nhit"] = cs1[cname][[1]] + cs2[cname][[1]]
      erchems[i,"nhit.agonist"] = cs1[cname][[1]]
      erchems[i,"nhit.antagonist"] = cs2[cname][[1]]
      vals = mat[,cname]
      vals = vals[vals>0]
      q = quantile(vals,probs=seq(0,1,0.1))
      if(length(vals)>0) erchems[i,"httr.pod"] = 3-q[10]
    }
  }

  temp = ersig0[is.element(ersig0$signature,sigs),]
  mat = reshape2::dcast(temp,signature~dtxsid,fill=0,value.var="logbmd",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  dmat = as.matrix(dist(t(mat)))
  ctemp = erchems[!is.na(erchems$auc.agonist),]
  d1 = ctemp[ctemp$auc.agonist>0.2,"dtxsid"]
  d2 = ctemp[ctemp$auc.antagonist>0.2,"dtxsid"]

  erchems$httr.mode = "not determined"
  for(i in 1:nrow(erchems)) {
    dtxsid = erchems[i,"dtxsid"]
    x1 = min(dmat[dtxsid,d1])
    x2 = min(dmat[dtxsid,d2])
    if(!is.na(x1) && !is.na(x2)) {
      if(x1<x2) erchems[i,"httr.mode"] = "agonist"
      else erchems[i,"httr.mode"] = "antagonist"
    }
    cat(erchems[i,"name"],erchems[i,"httr.mode"],"\n")
    #browser()
  }
  erchems[erchems$nhit<=5,"httr.mode"] = "inactive"

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
  mat = mat[,order(cs,decreasing=T)]
  names = colnames(mat)
  col.colors = names
  col.colors[] = "white"
  for(i in 1:length(names)) {
    name = names[i]
    if(is.element(name,supermat$newname)) {
      auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
      auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
      if(auc1>0.1 && auc2<0.1) col.colors[i] = "blue"
      else if(auc2>0.1 && auc1<0.1) col.colors[i] = "red"
      else if(auc2>0.1 && auc1>0.1) {
        if(auc1>auc2) col.colors[i] = "blue"
        else col.colors[i] = "red"
      }
      else col.colors[i] = "gray"
    }
  }

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,20),
                      dendrogram="row",
                      scale="none",
                      main="Reference Chemical\nSignatures Scaled Top",
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
                      ColSideColors=col.colors)
  chem.names = colnames(mat)
  sigorder = rev(rownames(mat)[result$rowInd])
  if(!to.file && do.stop) browser()

  ################################################################################################
  # BMD for the selected signatures and reference chemicals
  ################################################################################################
  temp = ersig0[is.element(ersig0$dtxsid,dlist),]
  chems = unique(temp[,c("dtxsid","name")])
  rownames(chems) = chems$dtxsid
  chems = chems[dlist,]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="logbmd",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat = mat[,chems$name]
  mat = mat[sigs,]
  mat = mat[,chem.names]
  mat = mat[sigorder,]

  names = colnames(mat)
  col.colors = names
  col.colors[] = "white"
  for(i in 1:length(names)) {
    name = names[i]
    if(is.element(name,supermat$newname)) {
      auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
      auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
      if(auc1>0.1 && auc2<0.1) col.colors[i] = "blue"
      else if(auc2>0.1 && auc1<0.1) col.colors[i] = "red"
      else if(auc2>0.1 && auc1>0.1) {
        if(auc1>auc2) col.colors[i] = "blue"
        else col.colors[i] = "red"
      }
      else col.colors[i] = "gray"
    }
  }

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,20),
                      dendrogram="none",
                      scale="none",
                      main="Reference Chemicals\nBMD Signatures Filtered",
                      xlab="",
                      ylab="",
                      cexCol=0.8,
                      cexRow=0.8,
                      Rowv=F,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="log(bmd)",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()

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
  mat = mat[,chem.names]
  limit = 10
  mat[mat> limit] = limit
  mat[mat< -limit] = -limit
  mat = mat[sigorder,]

  sig.names = rownames(mat)
  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,20),
                      dendrogram="none",
                      scale="none",
                      main="Reference Chemical\nSignatures AUC",
                      xlab="",
                      ylab="",
                      cexCol=0.8,
                      cexRow=0.8,
                      Rowv=F,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top directieon",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()

  ################################################################################################
  # hitcall for the selected signatures and reference chemicals
  ################################################################################################
  temp = ersig0[is.element(ersig0$dtxsid,dlist),]
  temp = temp[is.element(temp$signature,sigs),]
  chems = unique(temp[,c("dtxsid","name")])
  rownames(chems) = chems$dtxsid
  chems = chems[dlist,]

  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="hitcall",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat = mat[,chems$name]
  mat[is.na(mat)] = 0
  mat = mat[,chem.names]
  mat = mat[sigorder,]

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,20),
                      dendrogram="none",
                      scale="none",
                      main="Reference Chemical\nSignatures Hitcall",
                      xlab="",
                      ylab="",
                      cexCol=0.8,
                      cexRow=0.8,
                      Rowv=F,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="top directieon",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()

  ################################################################################################
  # select the chemicals that hit a minimum number of signatures and add to erchems$candidate
  ################################################################################################
  mask = erchems$candidate
  mask[] = 0
  mask[erchems$nhit.agonist>1] = 1
  mask[erchems$nhit.antagonist>1] = 1
  mask[erchems$httr.pod>2] = 0
  mask[is.element(erchems$dtxsid,refchems$dtxsid)] = 1
  y = erchems$auc.agonist
  y[is.na(y)] = 0
  mask[y>0.1] = 1
  y = erchems$auc.antagonist
  y[is.na(y)] = 0
  mask[y>0.1] = 1
  cat("selected chems 1:",sum(mask),"\n")
  erchems$candidate = mask
  candidates.dlist = erchems[erchems$candidate==1,"dtxsid"]
  candidates = erchems[erchems$candidate==1,"name"]
  file = paste0("../ERModel/ER_chems all ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  write.xlsx(erchems,file)

  ################################################################################################
  # candidates - scaled top
  ################################################################################################
  temp = ersig0[is.element(ersig0$signature,sigs),]
  temp = temp[is.element(temp$dtxsid,candidates.dlist),]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="scaled_top",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat[is.na(mat)] = 0
  mat = mat[,candidates]

  names = colnames(mat)
  col.colors = names
  col.colors[] = "white"
  for(i in 1:length(names)) {
    name = names[i]
    if(is.element(name,supermat$newname)) {
      auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
      auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
      if(auc1>0.1 && auc2<0.1) col.colors[i] = "blue"
      else if(auc2>0.1 && auc1<0.1) col.colors[i] = "red"
      else if(auc2>0.1 && auc1>0.1) {
        if(auc1>auc2) col.colors[i] = "blue"
        else col.colors[i] = "red"
      }
      else col.colors[i] = "gray"
    }
  }

  limit = 3
  mat[mat>limit] = limit
  mat[mat< -limit] = -limit

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,15),
                      dendrogram="row",
                      scale="none",
                      main="ER Candidates\nSignatures Scaled Top",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.8,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=col.colors)
  sigorder = rev(rownames(mat)[result$rowInd])

  ################################################################################################
  # Candidates - bmd
  ################################################################################################
  temp = ersig0[is.element(ersig0$signature,sigs),]
  temp = temp[is.element(temp$dtxsid,candidates.dlist),]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="logbmd",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,candidates]
  mat = mat[sigorder,]
  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,15),
                      dendrogram="none",
                      scale="none",
                      main="ER Candidates\nSignatures BMD",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.8,
                      Rowv=F,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()

  ################################################################################################
  # Get the ER actives based on the sig list - auc
  ################################################################################################
  temp = ersig0[is.element(ersig0$signature,sigs),]
  temp = temp[is.element(temp$dtxsid,candidates.dlist),]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="auc",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat[is.na(mat)] = 0
  mat = mat[,candidates]
  mat = mat[sigorder,]
  tmat = as.data.frame(t(mat))
  tmat = cbind(rownames(tmat),tmat,stringsAsFactors=F)
  tmat = cbind(rownames(tmat),tmat,stringsAsFactors=F)
  names(tmat)[1] = "dtxsid"
  names(tmat)[2] = "name"
  tmat = tmat[chems$name,]
  tmat$dtxsid = chems$dtxsid
  file = paste0("../ERModel/ER_chems_auc_signature ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  write.xlsx(tmat,file)

  limit = 10
  mat[mat>limit] = limit
  mat[mat< -limit] = -limit

  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,15),
                      dendrogram="none",
                      scale="none",
                      main="ER Candidates\nSignatures AUC",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.8,
                      Rowv=F,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()

  ################################################################################################
  # Candidates - hitcall
  ################################################################################################
  temp = ersig0[is.element(ersig0$signature,sigs),]
  temp = temp[is.element(temp$dtxsid,candidates.dlist),]
  mat = reshape2::dcast(temp,signature~name,fill=0,value.var="hitcall",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat[is.na(mat)] = 0
  mat = mat[,candidates]
  mat = mat[sigorder,]
  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,15),
                      dendrogram="none",
                      scale="none",
                      main="ER Candidates\nSignatures Hitcall",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.8,
                      Rowv=F,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()

  ################################################################################################
  # pull out the interesting genes - BMD
  ################################################################################################
  cat("pull out the interesting genes - BMD\n")
  genes = NULL
  for(i in 1:length(sigs)) {
    sig = sigs[i]
    temp = catalog[is.element(catalog$parent,sig),"signature"]
    for(x in temp) {
      y = genelists[x][[1]]
      genes = c(genes,y)
    }
  }
  allgenes = genes
  genetable = table(allgenes)
  genes = unique(genes)
  cat("Number of genes in signatures:",length(genes),"\n")
  genecr = genecr[is.element(genecr$dtxsid,candidates.dlist),]
  genecr = genecr[is.element(genecr$gene,genes),]
  genecr = genecr[genecr$top_over_cutoff>=0.5,]
  mat = reshape2::dcast(genecr,gene~name,fill=0,value.var="logbmd",fun.aggregate=mean)
  gcandidates = candidates[is.element(candidates,names(mat))]
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat = mat[,gcandidates]
  imat = abs(mat)
  imat[imat>0] = 1
  rs = rowSums(imat)
  mask = imat[,1]
  mask[] = 0
  mask[rs>=25] = 1
  mat = mat[mask==1,]
  selected.genes = rownames(mat)
  cat("Number of selected genes in signatures:",length(selected.genes),"\n")
  names = colnames(mat)
  col.colors = names
  col.colors[] = "white"
  for(i in 1:length(names)) {
    name = names[i]
    if(is.element(name,supermat$newname)) {
      auc1 = supermat[is.element(supermat$newname,name),"AUC.Agonist"][1]
      auc2 = supermat[is.element(supermat$newname,name),"AUC.Antagonist"][1]
      if(auc1>0.1 && auc2<0.1) col.colors[i] = "blue"
      else if(auc2>0.1 && auc1<0.1) col.colors[i] = "red"
      else if(auc2>0.1 && auc1>0.1) {
        if(auc1>auc2) col.colors[i] = "blue"
        else col.colors[i] = "red"
      }
      else col.colors[i] = "gray"
    }
  }
  print(dim(mat))
  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,15),
                      dendrogram="row",
                      scale="none",
                      main="ER Candidates\nGene BMD",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.2,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()

  ################################################################################################
  # pull out the interesting genes - scaled top
  ################################################################################################
  cat("pull out the interesting genes - scale top\n")
  genecr = genecr0
  genecr = genecr[is.element(genecr$dtxsid,candidates.dlist),]
  genecr = genecr[is.element(genecr$gene,selected.genes),]
  mat = reshape2::dcast(genecr,gene~name,fill=0,value.var="scaled_top",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat = mat[,gcandidates]
  mat[is.na(mat)] = 0
  imat = abs(mat)
  imat[imat>0] = 1
  rs = rowSums(imat)
  mask = imat[,1]
  mask[] = 0
  mask[rs>=10] = 1
  mat = mat[mask==1,]

  limit = 2
  mat[mat>limit] = limit
  mat[mat< -limit] = -limit
  print(dim(mat))
  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,15),
                      dendrogram="row",
                      scale="none",
                      main="ER Candidates\nGene Scaled Top",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.2,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()
  mat0 = mat
  ################################################################################################
  # classify genes as different between agonist and antagonist
  ################################################################################################
  genecr = genecr0
  genecr = genecr[is.element(genecr$dtxsid,candidates.dlist),]
  genecr = genecr[is.element(genecr$gene,selected.genes),]
  mat = reshape2::dcast(genecr,gene~dtxsid,fill=0,value.var="scaled_top",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  #dcandidates = chems[is.element(chems$name,candidates),"dtxsid"]
  mat = mat[,candidates.dlist]
  mat[is.na(mat)] = 0
  genelist = rownames(mat)
  name.list = c("gene","mean.agonist","mean.antagonist","diff")
  generes = as.data.frame(matrix(nrow=length(genelist),ncol=length(name.list)))
  names(generes) = name.list
  generes$gene = genelist
  rownames(generes) = generes$gene
  ctemp = erchems[!is.na(erchems$auc.agonist),]
  d1 = ctemp[ctemp$auc.agonist>0.4,"dtxsid"]
  d2 = ctemp[ctemp$auc.antagonist>0.4,"dtxsid"]
  d1 = d1[is.element(d1,names(mat))]
  d2 = d2[is.element(d2,names(mat))]
  generes$nsig = 0
  for(i in 1:nrow(generes)) {
    gene = generes[i,"gene"]
    x1 = mean(as.numeric(mat[gene,d1]))
    x2 = mean(as.numeric(mat[gene,d2]))
    generes[i,"mean.agonist"] = x1
    generes[i,"mean.antagonist"] = x2
    generes[i,"diff"] = x1-x2
    if(is.element(gene,names(genetable))) generes[i,"nsig"] = genetable[gene][[1]]
  }
#browser()
  x = generes$nsig
  y = generes$diff
  mod = lm(y~x)
  smod = summary(mod)
  intercept = smod$coefficients[1,1]
  slope = smod$coefficients[2,1]
  p = smod$coefficients[2,4]
  minsig = 5
  plot(y~x,xlab="Number of Signatures",ylab="Agonist-Antagonist Difference",
       ylim=c(-3,5),xlim=c(0,40),cex.lab=1.2,cex.axis=1.2,main="T/C difference vs. nsig")
  lines(c(0,50),c(intercept,intercept+slope*50))
  lines(c(0,100),c(0,0))
  lines(c(0,100),c(1,1))
  lines(c(0,100),c(-1,-1))
  lines(c(minsig,minsig),c(-10,10))
  text(30,-2,paste("p=",format(p,digits=2)),pos=4)

  if(!to.file && do.stop) browser()
  file = paste0("../ERModel/ER_gene_diff ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  write.xlsx(generes,file)
  mask = generes$nsig
  mask[] = 1
  mask[abs(generes$diff)<1] = 0
  mask[generes$nsig<minsig] = 0
  diffgenes = generes[mask==1,"gene"]
  mat2 = mat0[diffgenes,]
  result <- heatmap.2(as.matrix(mat2),
                      margins=c(10,15),
                      dendrogram="row",
                      scale="none",
                      main="ER Candidates\nSelected Gene Scaled Top",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.8,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()

  ################################################################################################
  # pull out the interesting genes - auc
  ################################################################################################
  cat("pull out the interesting genes - AUC\n")
  genecr = genecr0
  genecr = genecr[is.element(genecr$dtxsid,erchems$dtxsid),]
  genecr = genecr[is.element(genecr$gene,selected.genes),]
  mat = reshape2::dcast(genecr,gene~name,fill=0,value.var="auc",fun.aggregate=mean)
  rownames(mat) = mat[,1]
  mat = mat[,2:ncol(mat)]
  mat[is.na(mat)] = 0
  #browser()
  mat = mat[,gcandidates]
  imat = abs(mat)
  imat[imat>0] = 1
  rs = rowSums(imat)
  mask = imat[,1]
  mask[] = 0
  mask[rs>=10] = 1
  mat = mat[mask==1,]

  tmat = as.data.frame(t(mat))
  tmat = cbind(rownames(tmat),tmat,stringsAsFactors=F)
  tmat = cbind(rownames(tmat),tmat,stringsAsFactors=F)
  names(tmat)[1] = "dtxsid"
  names(tmat)[2] = "name"
  tmat = tmat[chems$name,]
  tmat$dtxsid = chems$dtxsid
  file = paste0("../ERModel/ER_chems_auc_gene ",dataset," ",sigset," ",hccut," ",minhit,".xlsx")
  write.xlsx(tmat,file)

  limit = 5
  mat[mat>limit] = limit
  mat[mat< -limit] = -limit
  print(dim(mat))
  result <- heatmap.2(as.matrix(mat),
                      margins=c(10,15),
                      dendrogram="row",
                      scale="none",
                      main="ER Candidates\nGene AUC",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.2,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()
  mat2 = mat[diffgenes,]
  result <- heatmap.2(as.matrix(mat2),
                      margins=c(10,15),
                      dendrogram="row",
                      scale="none",
                      main="ER Candidates\nSelected Gene AUC",
                      xlab="",
                      ylab="",
                      cexCol=0.3,
                      cexRow=0.8,
                      Rowv=T,
                      Colv=F,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(7,"RdGy"),
                      key.title="Key",
                      key.xlab="top direction",
                      cex.main=1,
                      ColSideColors=col.colors)
  if(!to.file && do.stop) browser()
  ################################################################################################
  # Plot ER model vs signature stats
  ################################################################################################
  cat("Plot ER model vs signature stats\n")
  par(mfrow=c(1,1),mar=c(4,4,4,2))
  temp = erchems[!is.na(erchems$auc.agonist),]
  temp = temp[!is.na(temp$auc.antagonist),]
  x = temp$auc.agonist
  for(i in 1:nrow(temp)) if(temp[i,"auc.antagonist"]>temp[i,"auc.agonist"]) x[i] = temp[i,"auc.antagonist"]
  y = temp$httr.pod

  mod = lm(y~x)
  smod = summary(mod)
  intercept = smod$coefficients[1,1]
  slope = smod$coefficients[2,1]
  p = smod$coefficients[2,4]
  plot(y~x,xlab="max(AUC(agonist),AUC(antagonist))",ylab="Mean(log(bmd uM))",cex.lab=1.2,cex.axis=1.2,ylim=c(-3,3),main=paste("p=",format(p,digits=2)))
  for(i in 1:nrow(temp)) {
    ag = temp[i,"auc.agonist"]
    antag = temp[i,"auc.antagonist"]
    col = "blue"
    if(antag>ag) col = "red"
    points(x[i],y[i],pch=21,bg=col)
  }
  lines(c(0,1.5),c(intercept,intercept+slope*1.5))

  if(!to.file) browser()
  else dev.off()
}

