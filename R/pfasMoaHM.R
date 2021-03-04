#--------------------------------------------------------------------------------------
#'
#' Build a heatmap of the stress genes
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
pfasMoaHM <- function(to.file=F,
                      sigset="screen_large",
                      method="fc",
                      hccut=0.95,
                      tccut=1.5) {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/PFAS/pfasMoaHM.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  dataset = "PFAS_HepaRG"
  celltype = "HepaRG"
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_active.xlsx")
  heparg = read.xlsx(file)

  dataset = "PFAS_U2OS"
  celltype = "U2OS"
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_active.xlsx")
  u2os = read.xlsx(file)

  file = "../input/signatures/super_target_key_for_pfas.xlsx"
  index = read.xlsx(file)

  file = "../input/pfas/PFAS categories.xlsx"
  cats = read.xlsx(file)
  cats$name = ""

  chems= unique(rbind(u2os[,c("dtxsid","name")],heparg[,c("dtxsid","name")]))
  chems$category = "-"
  for(i in 1:nrow(chems)) {
    dtxsid = chems[i,"dtxsid"]
    temp = sort(cats[is.element(cats$dtxsid,dtxsid),"category"])
    temp = paste(temp,collapse="|")
    chems[i,"category"] = temp
  }

  for(celltype in c("U2OS","HepaRG")) {
    for(tclass in unique(index$target_class)) {
      cat(tclass,"\n")
      st.list = index[is.element(index$target_class,tclass),"super_target"]
      if(celltype=="HepaRG") temp = heparg[is.element(heparg$super_target,st.list),]
      if(celltype=="U2OS") temp = u2os[is.element(u2os$super_target,st.list),]

      temp = temp[,c("name","super_target","bmd_median")]
      mat = reshape2::dcast(temp,name~super_target,value.var="bmd_median",fill=1000,fun.aggregate=median)
      rownames(mat) = mat[,1]
      mat = mat[,2:ncol(mat)]
      mat = 3-log10(mat)
      cex = 1
      if(ncol(mat)>20) cex=0.5
      result <- heatmap.2(as.matrix(mat),
                          margins=c(10,15),
                          scale="none",
                          main=paste(celltype,":",tclass),
                          xlab="",
                          ylab="",
                          cexCol=cex,
                          cexRow=0.5,
                          col=brewer.pal(9,"Reds"),
                          Rowv=T,
                          Colv=T,
                          dendrogram="both",
                          trace="none",
                          hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                          key=T,
                          key.title="Key",
                          key.xlab="log(bmd)")

      imat = mat
      imat[imat>0] = 1
      imat$sum = rowSums(imat)
      imat$name = rownames(imat)
      imat$category="-"
      for(k in 1:nrow(imat)) {
        name = rownames(imat)[k]
        imat[k,"category"] = chems[is.element(chems$name,name),"category"]
      }
      tclass1=tclass
      tclass1 = str_replace(tclass1,"\\/","_")
      file = paste0("../output/PFAS/pfasMoaHM_",celltype,"_",tclass1,".xlsx")
      write.xlsx(imat,file)
      if(!to.file) browser()
    }
  }
  if(to.file) dev.off()
}

