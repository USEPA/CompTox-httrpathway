#--------------------------------------------------------------------------------------
#' Generate the PFAS input data set
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_toxcast_pfas_pe1_normal
#--------------------------------------------------------------------------------------
pfasPlotCR <- function(to.file=F,
                       do.load=F,
                       dataset="PFAS_U2OS",
                       sigset="screen_large",
                       method="fc",
                       celltype="U2OS",
                       hccut=0.95,
                       tccut=2,
                       bmdcut=1) {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/signature_pfas/",celltype,"/pfas_crplots_",dataset,"_",sigset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }
  mat = MAT

  mat = mat[mat$hitcall>=hccut,]
  mat = mat[mat$top_over_cutoff>=tccut,]
  mat = mat[mat$bmd<bmdcut,]

  res = reshape2::dcast(mat,name~super_target,fill=0,value.var="auc",fun.aggregate=mean)
  rownames(res) = res[,1]
  res = as.matrix(res[,2:ncol(res)])
  res1 = res
  res1[res1>0] = 1
  cs =colSums(res1)
  #res = res[,cs>1]
  print(dim(res1))

  result <- heatmap.2(res,
                      margins=c(10,15),
                      dendrogram="both",
                      scale="none",
                      xlab="",
                      ylab="",
                      cexCol=0.6,
                      cexRow=0.6,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="AUC",
                      cex.main=1,
                      main=paste(celltype,"x target"))


  if(!to.file) browser()

  #################################################################################
  par(mfrow=c(3,2),mar=c(4,4,2,2))
  mat$proper_name = mat$name
  name.list = sort(unique(mat$name))
  for(name in name.list) {
    cat(name,"\n")
    temp1 = mat[is.element(mat$name,name),]
    st.list <- sort(unique(temp1$super_target))
    for(st in st.list) {
      cat("  ",st,"\n")
      temp2 <- temp1[is.element(temp1$super_target,st),]
      if(nrow(temp2)>0) {
        temp3 = temp2[order(temp2$top,decreasing=T),]
        for(i in 1:nrow(temp2)){
          signatureConcRespPlot(temp2[i,])
          if(!to.file) browser()
        }
      }
    }
  }
  if(!to.file) browser()
  else dev.off()

}

