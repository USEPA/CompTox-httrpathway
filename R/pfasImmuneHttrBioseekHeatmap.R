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
pfasImmuneHttrBioseekHeatmap <- function(to.file=F,
                                         do.load=F,
                                         tccut=2) {
  printCurrentFunction()

  if(do.load) {
    file = "../input/PFAS/PFAS wide hts.xlsx"
    cat(file,"\n")
    hts = read.xlsx(file)
    cat(nrow(hts),"\n")
    hts = hts[is.element(hts$source,"BSK"),]
    cat(nrow(hts),"\n")
    hts = hts[hts$hitcall>0.9,]
    cat(nrow(hts),"\n")
    HTS <<- hts

    celltype = "U2OS"
    file = paste0("../output/PFAS/pfasImmuneSummary_specific_",celltype,".xlsx")
    cat(file,"\n")
    u2os = read.xlsx(file)
    celltype = "HepaRG"

    file = paste0("../output/PFAS/pfasImmuneSummary_specific_",celltype,".xlsx")
    cat(file,"\n")
    heparg = read.xlsx(file)
    U2OS <<- u2os
    HEPARG <<- heparg
  }
  hts = HTS
  for(i in 1:nrow(hts)) {
    aname = hts[i,"hts_assay_name"]
    top = hts[i,"top"]
    if(contains(aname,"down")) hts[i,"top"] = -top
  }
  heparg = HEPARG
  u2os = U2OS


  file = "../input/PFAS/immune reference chemicals.xlsx"
  phenotype = "immunosuppression strong"
  map = read.xlsx(file)
  map = map[is.element(map$phenotype,phenotype),]

  heparg[is.element(heparg$dtxsid,map$dtxsid),"cclass"] = "Refchem"
  u2os[is.element(u2os$dtxsid,map$dtxsid),"cclass"] = "Refchem"

  heparg = heparg[is.element(heparg$cclass,c("PFAS","Refchem")),]
  u2os = u2os[is.element(u2os$cclass,c("PFAS","Refchem")),]

  heparg.cmat = unique(heparg[,c("name","dtxsid","cclass")])
  u2os.cmat = unique(u2os[,c("name","dtxsid","cclass")])
  heparg = heparg[heparg$top_over_cutoff>tccut,]
  u2os = u2os[u2os$top_over_cutoff>tccut,]

  names(hts)[4] = "sample_id"
  names(hts)[5] = "signature"
  name.list = c("sample_id","dtxsid","name","signature","top")

  heparg = heparg[,name.list]
  u2os = u2os[,name.list]
  hts = hts[,name.list]

  if(to.file) {
    fname <- paste0("../output/PFAS/pfasImmuneHttrBioseekHeatmap.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  for(celltype in c("U2OS","HepaRG")) {

    if(celltype == "U2OS") temp = u2os
    else temp = heparg
    hts0 = hts[is.element(hts$name,temp$name),]
    temp = rbind(temp,hts0)
    mat = reshape2::dcast(temp,name~signature,value.var="top",fun.aggregate=max,fill=0)
    rownames(mat) = mat[,1]
    mat = mat[,2:ncol(mat)]
    mat = t(mat)
    mat[mat>0] = 1
    mat[mat<0] = -1
    amat = abs(mat)
    rs = rowSums(amat)
    mat = mat[rs>=2,]
    cat(nrow(mat),"\n")
    cexrow = 0.1
    colors=rownames(mat)
    colors[] = "white"
    for(i in 1:nrow(mat)) {
      rname = rownames(mat)[i]
      if(contains(rname,"BSK")) {
        if(contains(rname,"down")) colors[i] = "red"
        else if(contains(rname,"up")) colors[i] = "green"
      }
    }

    ccolors = colnames(mat)
    for(i in 1:ncol(mat)) {
      cname = colnames(mat)[i]
      if(celltype=="HepaRG") {
        cclass = heparg.cmat[is.element(heparg.cmat$name,cname),"cclass"]
        if(cclass=="PFAS") ccolors[i] = "cyan"
        if(cclass=="Refchem") ccolors[i] = "red"
      }
      if(celltype=="U2OS") {
        cclass = u2os.cmat[is.element(u2os.cmat$name,cname),"cclass"]
        if(cclass=="PFAS") ccolors[i] = "cyan"
        if(cclass=="Refchem") ccolors[i] = "red"
      }
    }

    result <- heatmap.2(as.matrix(mat),
                        margins=c(10,10),
                        dendrogram="both",
                        scale="none",
                        main=celltype,
                        xlab="",
                        ylab="",
                        cexCol=0.2,
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
                        RowSideColors=colors,
                        ColSideColors=ccolors)
    if(!to.file) browser()

  }
  if(to.file) dev.off()
}
