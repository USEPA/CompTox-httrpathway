#--------------------------------------------------------------------------------------
#' Compile the summary statistics for the super targets
#' select chemcial / super targets pairs with good support
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#'
#--------------------------------------------------------------------------------------
superTargetStatsSummaryHeatmap <- function(to.file=F,
                                           do.load=F,
                                           dataset="heparg2d_toxcast_pfas_pe1_normal",
                                           sigset="screen_large",
                                           method="fc",
                                           celltype="HepaRG",
                                           hccut=0.95,
                                           tccut=1.5,
                                           target_class="Cancer") {
  printCurrentFunction()

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    mat = mat[mat$hitcall>=hccut,]
    mat = mat[mat$top_over_cutoff>=tccut,]
    SMAT <<- mat
  }
  smat = SMAT
  file = "../input/signatures/signatureDB_master_catalog 2021-03-05.xlsx"
  catalog = read.xlsx(file)
  catalog = catalog[catalog[,sigset]==1,]
  catalog = catalog[catalog[,"target_class"]==target_class,]
  st.list = sort(unique(catalog$super_target))
  smat = smat[is.element(smat$super_target,st.list),]

  file = paste0("../output/super_target_boxplot/superTargetStatsSummary.xlsx")
  mat = read.xlsx(file)
  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/",celltype,"/superTargetStatsSummaryHeatmap ",celltype," ",target_class,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  ctemp = unique(mat[,c("name","chem_super_target","use_class"),])
  refchems = NULL
  name.list = c("name","level")
  row = as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(row) = name.list
  if(target_class=="Cancer") {
    file = "../input/cancer/NTP ICE Cancer Data.xlsx"
    temp = read.xlsx(file)
    refchems = temp[,c("name","ilevel")]
    names(refchems) = name.list
  }
  else if(target_class=="Immune") {
    for(i in 1:nrow(ctemp)) {
      if(contains(ctemp[i,"name"],"Cyclosporin")) {
        row[1,"name"] = ctemp[i,"name"]
        row[1,"level"] = 4
        refchems = rbind(refchems,row)
      }
      if(contains(ctemp[i,"use_class"],"PAH")) {
        row[1,"name"] = ctemp[i,"name"]
        row[1,"level"] = 4
        refchems = rbind(refchems,row)
      }
    }
  }
  rownames(refchems) = refchems$name
  for(st in st.list) {
    cat(st,"\n")
    tmat = mat[is.element(mat$super_target,st),]
    tmat = tmat[is.element(tmat$celltype,celltype),]
    dlist = unique(tmat$dtxsid)

    if(length(dlist)>=2) {

      stmat = smat[is.element(smat$super_target,st),]
      stmat = stmat[is.element(stmat$dtxsid,dlist),]

      res = reshape2::dcast(stmat,name~signature,value.var="top",fill=0,fun.aggregate=median)
      rownames(res) = res[,1]
      res = res[,2:ncol(res)]
      res[res<0] = -1
      res[res>0] = 1

      nref = 0
      for(i in 1:nrow(res)) {
        if(is.element(rownames(res)[i],refchems$name)) nref = nref+1
      }
      #browser()
      if(nref>=2) {
        dmat = as.matrix(dist(res))
        rc = refchems[is.element(refchems$name,rownames(dmat)),"name"]
        dmat0 = dmat[rc,]
        useme = rc
        if(nrow(dmat0)>1) {
          for(i in 1:nrow(dmat0)) {
            dmat0 = dmat0[,order(dmat0[i,])]
            nc = min(10,ncol(dmat0))
            nlist = colnames(dmat0)[1:nc]
            useme = c(useme,nlist)
          }
          useme = unique(useme)
          res = res[useme,]
          colors=rownames(res)
          colors[] = "white"
          for(i in 1:nrow(res)) {
            rname = rownames(res)[i]
            if(is.element(rname,refchems$name)) {
              level = refchems[rname,"level"]
              if(level==0) colors[i] = "green"
              else if(level==1) colors[i] = "yellow"
              else if(level==2) colors[i] = "orange"
              else if(level==3) colors[i] = "red"
            }
          }

          result <- heatmap.2(as.matrix(res),
                              margins=c(10,15),
                              dendrogram="both",
                              scale="none",
                              main=paste(st,"\n",celltype),
                              xlab="",
                              ylab="",
                              cexCol=0.2,
                              cexRow=0.5,
                              Rowv=T,
                              Colv=T,
                              trace="none",
                              hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                              key=T,
                              col=brewer.pal(3,"PRGn"),
                              key.title="Key",
                              key.xlab="log(bmd)",
                              cex.main=1,
                              RowSideColors=colors)
          if(!to.file) browser()
        }
      }
    }
  }

  if(to.file) dev.off()
}

