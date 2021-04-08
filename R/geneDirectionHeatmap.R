#--------------------------------------------------------------------------------------
#' Make a heatmap of genes x signature x immune chemicals
#'
#'  heparg2d_toxcast_pfas_pe1_normal
#'  u2os_toxcast_pfas_pe1_normal
#'
#'  PFAS_HepaRG
#'  PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
geneDirectionHeatmap <- function(to.file=F,
                                 do.load=F,
                                 celltype="HepaRG",
                                 dataset="heparg2d_toxcast_pfas_pe1_normal",
                                 pval=0.05,
                                 nametag="_conthits",
                                 sigcatalog="signatureDB_master_catalog 2021-03-05",
                                 sigset="screen_large",
                                 phenotype="immunosuppression strong") {
  printCurrentFunction(paste(phenotype))
  if(to.file) {
    fname <- paste0("../output/PFAS/geneDirectionHeatmap",phenotype,"_",celltype,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  catalog = read.xlsx(file)
  catalog = catalog[catalog[,sigset]==1,]
  catalog = catalog[is.element(catalog$target_class,"Immune"),]

  file = "../output/PFAS/pfasImmuneRefchemHeatmap_immunosuppression strong_signature_signatures_HepaRG.xlsx"
  sigs = read.xlsx(file)
  catalog = catalog[is.element(catalog$parent,sigs[,1]),]
  signames = unique(catalog[,c("signature","parent","super_target")])
  if(do.load) {
    file <- paste0("../output/gene_conc_resp_summary/GENE_CR_",dataset,"_", pval, nametag ,".RData")
    print(file)
    load(file=file)
    GENE_CR <<- GENE_CR
  }

  file = paste0("../output/PFAS/pfasImmuneSummary_specific_",celltype,".xlsx")
  heparg = read.xlsx(file)

  file = "../input/PFAS/Immuntox chemical evidence.xlsx"
  chems = read.xlsx(file)
  file = "../input/signatures/signatureDB_genelists.RData"
  load(file=file)
  gene_cr = GENE_CR[is.element(GENE_CR$name,chems$chemical),]
  allgenes = unique(gene_cr$gene)

  signames = signames[order(signames$signature),]
  for(i in 1:nrow(signames)) {
    sig = signames[i,"signature"]
    st = catalog[is.element(catalog$signature,sig),"super_target"]
    genelist = genelists[sig][[1]]
    genelist = genelist[is.element(genelist,allgenes)]
    cat(sig,length(genelist),"\n")
    temp = gene_cr[is.element(gene_cr$gene,genelist),]

    tmat = reshape2::dcast(temp,name~gene,value.var="top",fun.aggregate=max,fill=0)
    tcmat = reshape2::dcast(temp,name~gene,value.var="top_over_cutoff",fun.aggregate=max,fill=0)
    hcmat = reshape2::dcast(temp,name~gene,value.var="hitcall",fun.aggregate=max,fill=0)
    rownames(tmat) = tcmat[,1]
    tmat = tmat[,2:ncol(tmat)]
    rownames(tcmat) = tcmat[,1]
    tcmat = tcmat[,2:ncol(tcmat)]
    rownames(hcmat) = hcmat[,1]
    hcmat = hcmat[,2:ncol(hcmat)]

    tmat[is.na(tmat)] = 0
    tcmat[is.na(tcmat)] = 0
    hcmat[is.na(hcmat)] = 0

    hcmat[hcmat<0.25] = 0
    hcmat[hcmat>0] = 1
    tmat[tmat>0] = 1
    tmat[tmat<0] = -1
    mat = hcmat * tmat * tcmat
    mat = cbind(mat,mat[,1])
    names(mat)[ncol(mat)] = "signature"
    mat[,"signature"] = 0
    for(j in 1:nrow(mat)) {
      cname = rownames(mat)[j]
      temp = heparg[is.element(heparg$name,cname),]
      parent = signames[is.element(signames$signature,sig),"parent"]
      temp = temp[is.element(temp$signature,parent),]
      tc = 0
      if(nrow(temp) >0) {
        top = temp[1,"top"]
        tc =  temp[1,"top_over_cutoff"]
        hc =  temp[1,"hitcall"]
        if(top<0) tc = -1*tc
      }
      mat[j,"signature"] = tc
    }
    cutoff = 5
    mat[mat>cutoff] = cutoff
    mat[mat< -cutoff] = -cutoff
    result <- heatmap.2(as.matrix(mat),
                        margins=c(10,10),
                        dendrogram="both",
                        scale="none",
                        main=paste(sig,"\n",st),
                        xlab="",
                        ylab="",
                        cexCol=0.2,
                        cexRow=0.2,
                        Rowv=T,
                        Colv=T,
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        col=brewer.pal(9,"PRGn"),
                        key.title="Key",
                        key.xlab="T/C",
                        cex.main=0.5)
    if(!to.file) browser()
  }

  if(to.file) dev.off()

}
