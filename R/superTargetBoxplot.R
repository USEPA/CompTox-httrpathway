#--------------------------------------------------------------------------------------
#' Generate chemicalwise boxplot of the BMD distributions by super_target
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#'  u2os_toxcast_pfas_pe1_normal

#'
#--------------------------------------------------------------------------------------
superTargetBoxplot <- function(to.file=F,
                               do.load=F,
                               dataset="u2os_toxcast_pfas_pe1_normal",
                               sigset="screen_large",
                               method="fc",
                               celltype="U2OS",
                               hccut=0.95,
                               tccut=2.5) {
  printCurrentFunction(paste(dataset,sigset,method))

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

  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,15,2,2))
  chem.list = sort(unique(mat$name))
  for(chem in chem.list) {
    temp = mat[is.element(mat$name,chem),]
    cat(chem,nrow(temp),"\n")
    x = temp$super_target
    y = temp$bmd

    st.list = unique(x)

    res = as.data.frame(matrix(nrow=length(st.list),ncol=2))
    names(res) = c("super_target","bmd_median")
    res$super_target = st.list
    for(i in 1:length(st.list)) {
      st = st.list[i]
      temp = y[is.element(x,st)]
      res[i,2] = median(temp)
    }
    res$newname = ""
    res = res[order(res$bmd_median),]
    for(i in 1:nrow(res)) {
      st = res[i,"super_target"]
      si = as.character(i)
      if(i<100) si = paste0("0",i)
      if(i<10) si = paste0("00",i)
      res[i,"newname"] = paste(si, res[i,"super_target"])
    }
    xnew = x
    for(i in 1:nrow(res)) {
      oldname = res[i,"super_target"]
      newname = res[i,"newname"]
      xnew[is.element(x,oldname)] = newname
    }

    res$count = 0
    for(i in 1:nrow(res)) {
      newname = res[i,"newname"]
      res[i,"count"] = length(xnew[xnew==newname])
    }
    mask = y
    mask[] = 1
    exclude.list = res[res$count==1,"newname"]
    mask[is.element(xnew,exclude.list)] = 0
    #browser()

    xnew = xnew[mask==1]
    y = y[mask==1]

    res = res[is.element(res$newname,xnew),]
    nmax = 40
    if(nrow(res)>nmax) {
      res = res[1:nmax,]
      mask = xnew
      mask[] = 0
      mask[is.element(xnew,res$newname)] = 1
      xnew = xnew[mask==1]
      y = y[mask==1]
    }
    if(length(y)>0 && nrow(res)>1) {
      boxplot(y~xnew,main=paste(chem,"\n",celltype),
              ylim=c(0.01,100),log="x",xlab="BMD",ylab="",
              horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),names=res$super_target)
      for(bmd in c(100,10,1,0.1,0.001)) lines(c(bmd,bmd),c(0,1000000))
      if(!to.file) browser()
    }
  }
  if(to.file) dev.off()

}

