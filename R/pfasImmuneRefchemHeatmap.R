#--------------------------------------------------------------------------------------
#' Make a heatmap of the response from a reference chemica lset
#'
#'  heparg2d_toxcast_pfas_pe1_normal
#'  u2os_toxcast_pfas_pe1_normal
#'
#'  PFAS_HepaRG
#'  PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
pfasImmuneRefchemHeatmap <- function(to.file=F,
                                     do.load=F,
                                     sigcatalog="signatureDB_master_catalog 2021-03-05",
                                     sigset="screen_large",
                                     phenotype="immunosuppression strong",
                                     sclass="signature") {
  printCurrentFunction(paste(phenotype))

  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  catalog = read.xlsx(file)
  catalog = catalog[catalog[,sigset]==1,]
  signames = unique(catalog[,c("parent","super_target")])

  file = "../input/PFAS/immune reference chemicals.xlsx"
  map = read.xlsx(file)
  map = map[is.element(map$phenotype,phenotype),]
  celltype = "U2OS"
  file = paste0("../output/PFAS/pfasImmuneSummary_specific_",celltype,".xlsx")
  u2os = read.xlsx(file)
  celltype = "HepaRG"
  file = paste0("../output/PFAS/pfasImmuneSummary_specific_",celltype,".xlsx")
  heparg = read.xlsx(file)

  heparg = heparg[heparg$top_over_cutoff>2,]
  u2os = u2os[u2os$top_over_cutoff>2,]
  u2os0 = u2os
  heparg0 = heparg

  dlist = map$dtxsid
  u2os = u2os[is.element(u2os$dtxsid,dlist),]
  heparg = heparg[is.element(heparg$dtxsid,dlist),]
  target.names = unique(c(heparg$name,u2os$name))
  if(to.file) {
    fname <- paste0("../output/PFAS/pfasImmuneRefchemHeatmap_",phenotype,"_",sclass,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  for(celltype in c("U2OS","HepaRG")) {
  #for(celltype in c("HepaRG")) {

    if(celltype == "U2OS") {
      temp = u2os
      temp0 = u2os0
    }
    else {
      temp = heparg
      temp0 = heparg0
    }
    temp = temp[order(temp$super_target),]
    pfas.names = sort(unique(temp0[is.element(temp0$cclass,"PFAS"),"name"]))

    if(sclass=="super_target") mat = reshape2::dcast(temp,name~super_target,value.var="top",fun.aggregate=max,fill=0)
    if(sclass=="signature") mat = reshape2::dcast(temp,name~signature,value.var="top",fun.aggregate=max,fill=0)
    rownames(mat) = mat[,1]
    mat = mat[,2:ncol(mat)]
    mat = t(mat)
    mat[mat>0] = 1
    mat[mat<0] = -1
    amat = abs(mat)
    rs = rowSums(amat)
    mat = mat[rs>=2,]
    cat(nrow(mat),"\n")
    cexrow = 0.4
    if(sclass=="signature") cexrow = 0.1

    sigs = NULL
    rcolors = rownames(mat)
    for(i in 1:length(rcolors)) {
      sig = rcolors[i]
      st = signames[is.element(signames$parent,sig),"super_target"][1]
      rcolors[i] = st
    }

    stlist = rcolors
    rcolors[] = "white"
    for(i in 1:length(rcolors)) {
      if(stlist[i]=="B Cell") rcolors[i] = "red"
      if(stlist[i]=="B Cell vs Plasma Cell") rcolors[i] = "red"
      if(stlist[i]=="B Lymphocyte") rcolors[i] = "red"
      if(stlist[i]=="Dendritic Cell") rcolors[i] = "lightgray"
      if(stlist[i]=="Dendritic Cell vs Macrophage") rcolors[i] = "lightgray"
      if(stlist[i]=="Dendritic Cell vs Pathogen") rcolors[i] = "lightgray"

      if(stlist[i]=="Cytokine") rcolors[i] = "black"
      if(stlist[i]=="EGF") rcolors[i] = "black"
      if(stlist[i]=="EGF/EGFR") rcolors[i] = "black"
      if(stlist[i]=="NFKB") rcolors[i] = "black"
      if(stlist[i]=="TGFB") rcolors[i] = "black"

      if(stlist[i]=="Natural Killer Cell") rcolors[i] = "orange"
      if(stlist[i]=="Natural Killer Cell vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="Basophil vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="Memory T Cell") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell vs Thymocyte") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell") rcolors[i] = "cyan"
    }

    result <- heatmap.2(as.matrix(mat),
                        margins=c(10,10),
                        dendrogram="both",
                        scale="none",
                        main=paste(celltype,"\nrefchems"),
                        xlab="",
                        ylab="",
                        cexCol=1,
                        cexRow=cexrow,
                        Rowv=T,
                        Colv=T,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(3,"PRGn"),
                        key.title="Key",
                        key.xlab="top directieon",
                        cex.main=1,
                        RowSideColors=rcolors)
    if(!to.file) browser()
############################################################################
    sig.list = rownames(mat)
    if(celltype == "U2OS") temp = u2os0
    else temp = heparg0
    temp = temp[is.element(temp$signature,sig.list),]
    if(sclass=="super_target") mat = reshape2::dcast(temp,name~super_target,value.var="top",fun.aggregate=max,fill=0)
    if(sclass=="signature") mat = reshape2::dcast(temp,name~signature,value.var="top",fun.aggregate=max,fill=0)
    rownames(mat) = mat[,1]
    mat = mat[,2:ncol(mat)]
    mat = t(mat)
    mat[mat>0] = 1
    mat[mat<0] = -1
    amat = abs(mat)
    #rs = rowSums(amat)
    #mat = mat[rs>=2,]
    cs = colSums(amat)
    cutoff = 20
    if(celltype=="U2OS") cutoff = 10
    mat = mat[,cs>=cutoff]
    cat(dim(mat),"\n")
    cexrow = 0.1
    if(sclass=="signature") cexrow = 0.1

    chemnames = colnames(mat)
    colors = chemnames
    colors[] = "white"
    colors[is.element(chemnames,target.names)] = "red"
    colors[is.element(chemnames,pfas.names)] = "cyan"
    rcolors = rownames(mat)
    for(i in 1:length(rcolors)) {
      sig = rcolors[i]
      st = signames[is.element(signames$parent,sig),"super_target"][1]
      rcolors[i] = st
    }
    stlist = rcolors
    rcolors[] = "white"
    for(i in 1:length(rcolors)) {
      if(stlist[i]=="B Cell") rcolors[i] = "red"
      if(stlist[i]=="B Cell vs Plasma Cell") rcolors[i] = "red"
      if(stlist[i]=="B Lymphocyte") rcolors[i] = "red"
      if(stlist[i]=="Dendritic Cell") rcolors[i] = "lightgray"
      if(stlist[i]=="Dendritic Cell vs Macrophage") rcolors[i] = "lightgray"
      if(stlist[i]=="Dendritic Cell vs Pathogen") rcolors[i] = "lightgray"

      if(stlist[i]=="Cytokine") rcolors[i] = "black"
      if(stlist[i]=="EGF") rcolors[i] = "black"
      if(stlist[i]=="EGF/EGFR") rcolors[i] = "black"
      if(stlist[i]=="NFKB") rcolors[i] = "black"
      if(stlist[i]=="TGFB") rcolors[i] = "black"

      if(stlist[i]=="Natural Killer Cell") rcolors[i] = "orange"
      if(stlist[i]=="Natural Killer Cell vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="Basophil vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="Memory T Cell") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell vs Thymocyte") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell") rcolors[i] = "cyan"
    }
    result <- heatmap.2(as.matrix(mat),
                        margins=c(10,10),
                        dendrogram="both",
                        scale="none",
                        main=paste(celltype,"\nAll Chems"),
                        xlab="",
                        ylab="",
                        cexCol=0.1,
                        cexRow=cexrow,
                        Rowv=T,
                        Colv=T,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(3,"PRGn"),
                        key.title="Key",
                        key.xlab="top direction",
                        cex.main=1,
                        ColSideColors=colors,
                        RowSideColors=rcolors)
    if(!to.file) browser()

    ############################################################################
    sig.list = rownames(mat)
    if(celltype == "U2OS") temp = u2os0
    else temp = heparg0
    temp = temp[is.element(temp$signature,sig.list),]
    if(sclass=="super_target") mat = reshape2::dcast(temp,name~super_target,value.var="top",fun.aggregate=max,fill=0)
    if(sclass=="signature") mat = reshape2::dcast(temp,name~signature,value.var="top",fun.aggregate=max,fill=0)
    rownames(mat) = mat[,1]
    mat = mat[,2:ncol(mat)]
    mat = t(mat)
    mat[mat>0] = 1
    mat[mat<0] = -1
    amat = abs(mat)
    cs = colSums(amat)
    cutoff = 50
    if(celltype=="U2OS") cutoff = 50
    mask = cs
    mask[] = 0
    mask[cs>=cutoff] = 1
    cn = colnames(mat)
    mask[is.element(cn,pfas.names)] = 1
    mask[is.element(cn,target.names)] = 1

    mat = mat[,mask==1]
    cat(dim(mat),"\n")
    cexrow = 0.1

    chemnames = colnames(mat)
    colors = chemnames
    colors[] = "white"
    colors[is.element(chemnames,target.names)] = "red"
    colors[is.element(chemnames,pfas.names)] = "cyan"

    rcolors = rownames(mat)
    for(i in 1:length(rcolors)) {
      sig = rcolors[i]
      st = signames[is.element(signames$parent,sig),"super_target"][1]
      rcolors[i] = st
    }
    stlist = rcolors
    rcolors[] = "white"
    for(i in 1:length(rcolors)) {
      if(stlist[i]=="B Cell") rcolors[i] = "red"
      if(stlist[i]=="B Cell vs Plasma Cell") rcolors[i] = "red"
      if(stlist[i]=="B Lymphocyte") rcolors[i] = "red"
      if(stlist[i]=="Dendritic Cell") rcolors[i] = "lightgray"
      if(stlist[i]=="Dendritic Cell vs Macrophage") rcolors[i] = "lightgray"
      if(stlist[i]=="Dendritic Cell vs Pathogen") rcolors[i] = "lightgray"

      if(stlist[i]=="Cytokine") rcolors[i] = "black"
      if(stlist[i]=="EGF") rcolors[i] = "black"
      if(stlist[i]=="EGF/EGFR") rcolors[i] = "black"
      if(stlist[i]=="NFKB") rcolors[i] = "black"
      if(stlist[i]=="TGFB") rcolors[i] = "black"

      if(stlist[i]=="Natural Killer Cell") rcolors[i] = "orange"
      if(stlist[i]=="Natural Killer Cell vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="Basophil vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="Memory T Cell") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell vs Thymocyte") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell") rcolors[i] = "cyan"
    }
    result <- heatmap.2(as.matrix(mat),
                        margins=c(15,10),
                        dendrogram="both",
                        scale="none",
                        main=paste(celltype,"All Chems\n50 sigs+PFAS+Refchems"),
                        xlab="",
                        ylab="",
                        cexCol=0.1,
                        cexRow=cexrow,
                        Rowv=T,
                        Colv=T,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(3,"PRGn"),
                        key.title="Key",
                        key.xlab="top direction",
                        cex.main=1,
                        ColSideColors=colors,
                        RowSideColors=rcolors)
    temp = as.data.frame(matrix(nrow=nrow(mat),ncol=2))
    names(temp) = c("signature","celltype")
    temp[,1] = rownames(mat)
    temp[,2] = celltype
    sigs = rbind(sigs,temp)
    if(!to.file) browser()

    ############################################################################
    sig.list = rownames(mat)
    if(celltype == "U2OS") temp = u2os0
    else temp = heparg0
    temp = temp[is.element(temp$signature,sig.list),]
    if(sclass=="super_target") mat = reshape2::dcast(temp,name~super_target,value.var="top",fun.aggregate=max,fill=0)
    if(sclass=="signature") mat = reshape2::dcast(temp,name~signature,value.var="top",fun.aggregate=max,fill=0)
    rownames(mat) = mat[,1]
    mat = mat[,2:ncol(mat)]
    mat = t(mat)
    mat[mat>0] = 1
    mat[mat<0] = -1
    amat = abs(mat)
    cs = colSums(amat)
    cutoff = 50
    if(celltype=="U2OS") cutoff = 50
    mask = cs
    mask[] = 0
    #mask[cs>=cutoff] = 1
    cn = colnames(mat)
    mask[is.element(cn,pfas.names)] = 1
    mask[is.element(cn,target.names)] = 1

    mat = mat[,mask==1]
    cat(dim(mat),"\n")
    cexrow = 0.1

    chemnames = colnames(mat)
    colors = chemnames
    colors[] = "white"
    colors[is.element(chemnames,target.names)] = "red"
    colors[is.element(chemnames,pfas.names)] = "cyan"
    rcolors = rownames(mat)
    for(i in 1:length(rcolors)) {
      sig = rcolors[i]
      st = signames[is.element(signames$parent,sig),"super_target"][1]
      rcolors[i] = st
    }
    stlist = rcolors
    rcolors[] = "white"
    for(i in 1:length(rcolors)) {
      if(stlist[i]=="B Cell") rcolors[i] = "red"
      if(stlist[i]=="B Cell vs Plasma Cell") rcolors[i] = "red"
      if(stlist[i]=="B Lymphocyte") rcolors[i] = "red"
      if(stlist[i]=="Dendritic Cell") rcolors[i] = "lightgray"
      if(stlist[i]=="Dendritic Cell vs Macrophage") rcolors[i] = "lightgray"
      if(stlist[i]=="Dendritic Cell vs Pathogen") rcolors[i] = "lightgray"

      if(stlist[i]=="Cytokine") rcolors[i] = "black"
      if(stlist[i]=="EGF") rcolors[i] = "black"
      if(stlist[i]=="EGF/EGFR") rcolors[i] = "black"
      if(stlist[i]=="NFKB") rcolors[i] = "black"
      if(stlist[i]=="TGFB") rcolors[i] = "black"

      if(stlist[i]=="Natural Killer Cell") rcolors[i] = "orange"
      if(stlist[i]=="Natural Killer Cell vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="Basophil vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="Memory T Cell") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell vs Th Cell") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell vs Thymocyte") rcolors[i] = "cyan"
      if(stlist[i]=="T Cell") rcolors[i] = "cyan"
    }
    result <- heatmap.2(as.matrix(mat),
                        margins=c(15,10),
                        dendrogram="both",
                        scale="none",
                        main=paste(celltype,"All Chems\nPFAS+Refchems"),
                        xlab="",
                        ylab="",
                        cexCol=0.1,
                        cexRow=cexrow,
                        Rowv=T,
                        Colv=T,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(3,"PRGn"),
                        key.title="Key",
                        key.xlab="top direction",
                        cex.main=1,
                        ColSideColors=colors,
                        RowSideColors=rcolors)
    if(!to.file) browser()
    file <- paste0("../output/PFAS/pfasImmuneRefchemHeatmap_",phenotype,"_",sclass,"_signatures_",celltype,".xlsx")
    write.xlsx(sigs,file)
  }
  if(to.file) dev.off()

}
