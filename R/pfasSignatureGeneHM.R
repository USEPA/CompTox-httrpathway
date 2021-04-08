#--------------------------------------------------------------------------------------
#'
#' build a heatmap of the signatures showing string PFAS activity in immune
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
pfasSignatureGeneHM <- function(to.file=F,
                                do.load=F,
                                sigcatalog="signatureDB_master_catalog 2021-03-05",
                                sigset="screen_large",
                                phenotype="immunosuppression strong",
                                dataset="DMEM_6hr_pilot_normal_pe_1") {
  printCurrentFunction()

  if(to.file) {
    file = paste0("../output/PFAS/pfasImmuneRefchemHeatmap_",phenotype,"_signature_HM.pdf")
    pdf(file=file,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  if(do.load) {
    file = "../input/fcdata/FCMAT2_heparg2d_toxcast_pfas_pe1_normal.RData"
    load(file=file)
    heparg.xgenes <<- colnames(FCMAT2)

    file = "../input/fcdata/FCMAT2_u2os_toxcast_pfas_pe1_normal.RData"
    load(file=file)
    u2os.xgenes <<- colnames(FCMAT2)
  }

  file = paste0("../output/PFAS/pfasImmuneRefchemHeatmap_",phenotype,"_signature_signatures.xlsx")

  file = "../input/signatures/signatureDB_genelists.RData"
  load(file=file)

  file=paste0("../input/signatures/",sigcatalog,".xlsx")
  catalog = read.xlsx(file)
  catalog = catalog[catalog[,sigset]==1,]

  for(celltype in c("U2OS","HepaRG")) {
    file <- paste0("../output/PFAS/pfasImmuneRefchemHeatmap_",phenotype,"_signature_signatures_",celltype,".xlsx")
    sigs = read.xlsx(file)
    if(celltype=="U2OS") {
      sigs = unique(sigs[,"signature"])
      sigs = catalog[is.element(catalog$parent,sigs),"signature"]
      xgenes = u2os.genes
    }
    if(celltype=="HepaRG") {
      sigs = unique(sigs[,"signature"])
      sigs = catalog[is.element(catalog$parent,sigs),"signature"]
      xgenes = heparg.xgenes
     }

    genes = NULL

    for(i in 1:length(sigs)) {
      sig = sigs[i]
      genes = c(genes,genelists[sig][[1]])
      #cat(sig,length(genelists[sig][[1]]),"\n")
    }
    genes = unique(genes)
    genes = genes[is.element(genes,xgenes)]

    mat = as.data.frame(matrix(nrow=length(genes),ncol=length(sigs)))
    rownames(mat) = genes
    colnames(mat) = sigs
    mat[] = 0

    for(i in 1:length(sigs)) {
      sig = sigs[i]
      gl = genelists[sig][[1]]
      gl = gl[is.element(gl,xgenes)]
      mat[gl,sig] = 1
    }

    file = paste0("../output/PFAS/pfasImmuneRefchemHeatmap_",phenotype,"_signature_genecounts_",celltype,".xlsx")
    rs = rowSums(mat)
    nmax = 5
    if(celltype=="HepaRG") nmax = 20
    mat = mat[rs>=nmax,]
    x = as.data.frame(rs)
    x[,2] = rownames(x)
    x = x[,c(2,1)]
    names(x) = c("gene","count")
    write.xlsx(x,file)
    cs = colSums(mat)
    mat = mat[,cs>=20]

    colors = colnames(mat)
    for(i in 1:length(colors)) {
      sig = colors[i]
      st = catalog[is.element(catalog$signature,sig),"super_target"][1]
      colors[i] = st
    }
    stlist = colors
    colors[] = "white"
    for(i in 1:length(colors)) {
      if(stlist[i]=="B Cell") colors[i] = "red"
      else if(stlist[i]=="B Cell vs Plasma Cell") colors[i] = "red"
      else if(stlist[i]=="B Lymphocyte") colors[i] = "red"
      else if(stlist[i]=="Lymphocyte") colors[i] = "red"
      else if(stlist[i]=="B Cell vs Myeloid Cell") colors[i] = "red"
      else if(stlist[i]=="B Lymphocyte vs Plasma Cell") colors[i] = "red"

      else if(stlist[i]=="Dendritic Cell") colors[i] = "lightgray"
      else if(stlist[i]=="Dendritic Cell vs Macrophage") colors[i] = "lightgray"
      else if(stlist[i]=="Dendritic Cell vs Pathogen") colors[i] = "lightgray"
      else if(stlist[i]=="Dendritic Cell vs Monocyte") colors[i] = "lightgray"
      else if(stlist[i]=="Myeloid Dendritic Cell") colors[i] = "lightgray"


      else if(stlist[i]=="Cytokine") colors[i] = "black"
      else if(stlist[i]=="EGF") colors[i] = "black"
      else if(stlist[i]=="EGF/EGFR") colors[i] = "black"
      else if(stlist[i]=="NFKB") colors[i] = "black"
      else if(stlist[i]=="TGFB") colors[i] = "black"
      else if(stlist[i]=="TNF") colors[i] = "black"
      else if(stlist[i]=="Interferon") colors[i] = "black"
      else if(stlist[i]=="Interleukin") colors[i] = "black"
      else if(stlist[i]=="IFNG") colors[i] = "black"


      else if(stlist[i]=="Natural Killer Cell") colors[i] = "orange"

      else if(stlist[i]=="Natural Killer Cell vs Th Cell") colors[i] = "cyan"
      else if(stlist[i]=="Basophil vs Th Cell") colors[i] = "cyan"
      else if(stlist[i]=="Memory T Cell") colors[i] = "cyan"
      else if(stlist[i]=="T Cell vs Th Cell") colors[i] = "cyan"
      else if(stlist[i]=="T Cell vs Thymocyte") colors[i] = "cyan"
      else if(stlist[i]=="T Cell") colors[i] = "cyan"
      else if(stlist[i]=="Monocyte vs T Cell ") colors[i] = "cyan"
      else if(stlist[i]=="T Lymphocyte") colors[i] = "cyan"
      else if(stlist[i]=="Natural Killer Cell vs T Cell") colors[i] = "cyan"
      else if(stlist[i]=="Mast Cell vs Th Cell") colors[i] = "cyan"
      else if(stlist[i]=="B Cell vs T Cell") colors[i] = "cyan"
      else if(stlist[i]=="B Cell vs Th Cell") colors[i] = "cyan"
      else if(stlist[i]=="Dendritic Cell vs Th Cell") colors[i] = "cyan"
      else if(stlist[i]=="Lymphocyte vs T Cell") colors[i] = "cyan"
      else if(stlist[i]=="Macrophage vs Th Cell") colors[i] = "cyan"
      else if(stlist[i]=="Monocyte vs T Cell") colors[i] = "cyan"

      else if(stlist[i]=="PBMC") colors[i] = "khaki"
      else if(stlist[i]=="PBMC vs Pathogen") colors[i] = "khaki"
      else if(stlist[i]=="PBMC vs Vaccine") colors[i] = "khaki"
      else if(stlist[i]=="Monocyte vs Pathogen") colors[i] = "khaki"
      else if(stlist[i]=="Monocyte") colors[i] = "khaki"

      else if(stlist[i]=="Macrophage vs Interferon") colors[i] = "blue"
      else if(stlist[i]=="Macrophage vs Pathogen") colors[i] = "blue"
      else if(stlist[i]=="Macrophage vs Monocyte") colors[i] = "blue"
      else if(stlist[i]=="Macrophage") colors[i] = "blue"
      else if(stlist[i]=="Macrophage vs Interleukin") colors[i] = "blue"

      else if(stlist[i]=="Endothelial Cell vs Interferon") colors[i] = "yellow"
      else if(stlist[i]=="Complement") colors[i] = "yellow"
      else if(stlist[i]=="Fibroblast vs Interferon") colors[i] = "yellow"
      else if(stlist[i]=="Virus") colors[i] = "yellow"
      else if(stlist[i]=="Immune Disease") colors[i] = "yellow"
      else if(stlist[i]=="Thymocyte") colors[i] = "yellow"
      else if(stlist[i]=="Leukocyte vs Pathogen") colors[i] = "yellow"
      else if(stlist[i]=="Microglia Interferon") colors[i] = "yellow"
      else if(stlist[i]=="Arthritis") colors[i] = "yellow"
      else if(stlist[i]=="Chronic Lymphocytic Leukemia") colors[i] = "yellow"
      else if(stlist[i]=="Fibroblast vs Pathogen") colors[i] = "yellow"
      else if(stlist[i]=="Inflammation") colors[i] = "yellow"
      else if(stlist[i]=="Lupus") colors[i] = "yellow"
      else if(stlist[i]=="Mast Cell Leukemia") colors[i] = "yellow"


      else cat(stlist[i],"\n")
    }

    result <- heatmap.2(as.matrix(mat),
                        margins=c(5,5),
                        scale="none",
                        main=paste("PFAS Immune Signatures\n",celltype),
                        xlab="",
                        ylab="",
                        cexCol=0.1,
                        cexRow=0.3,
                        col=brewer.pal(3,"Reds"),
                        Rowv=T,
                        Colv=T,
                        dendrogram="both",
                        trace="none",
                        hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                        key=T,
                        key.title="Key",
                        key.xlab="l2fc",
                        cex.main=1,
                        ColSideColors=colors)

    if(!to.file) browser()
  }
  if(to.file) dev.off()
}

