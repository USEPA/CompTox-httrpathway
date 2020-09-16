#--------------------------------------------------------------------------------------
#' Generate chemicalwise boxplot of the BMD distributions by super_target
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#'  u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
superTargetBoxplot <- function(to.file=F,
                               do.load=F,
                               dataset="PFAS_HepaRG",
                               sigset="screen_large",
                               method="fc",
                               celltype="HepaRG",
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

  file = "../input/signatures/refchemdb_chem_filtered_unique_with_new.xlsx"
  file = "../input/chemicals/httr_chemical_annotations 2020-09-14.xlsx"
  rcdb = read.xlsx(file)
  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,15,6,2))

  name.list = c("name","celltype","dtxsid","sid","use","annotation","pod","burst_pod","ratio","specific_targets")
  row = as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(row) = name.list
  summary = NULL

  chem.list = sort(unique(mat$name))
  #chem.list = "Zearalenone"
  #chem.list = chem.list[1:10]
  for(chem in chem.list) {
    temp0 = mat[is.element(mat$name,chem),]
    dtxsid = temp0[1,"dtxsid"]
    chem_annotation = NULL
    chem_use = NULL
    if(is.element(dtxsid,rcdb$dtxsid)) {
      chem_annotation = unique(rcdb[is.element(rcdb$dtxsid,dtxsid),"annotation"])
      chem_use = unique(rcdb[is.element(rcdb$dtxsid,dtxsid),"use_class"])
    }
    sid.list = unique(temp0$sample_id)
    for(sid in sid.list) {
      temp = temp0[is.element(temp0$sample_id,sid),]

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
      burst_pod = 1000
      pod = 1000
      targets = "None"
      pod_target = "None"
      if(nrow(res)>1) {
        burst_pod = 10**(median(log10(res$bmd_median)))
        pod = min(res$bmd_median)
        pod_target = res[1,"super_target"]
        temp = res[res$bmd_median<burst_pod/10,]
        targets = "None"
        if(nrow(temp)>0) targets = paste(temp$super_target,collapse="|")
      }

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
        newnames = res$super_target
        cols = newnames
        cols[] = "white"
        for(k in 1:length(newnames)) {
          nnk = newnames[k]
          if(is.element(nnk,chem_annotation)) {
            nnk = paste(nnk,"*")
            cols[k] = "red"
          }
          newnames[k] = nnk
        }
        suse = paste(chem_use,collapse="|")
        star = paste(chem_annotation,collapse="|")
        main=paste(chem,"\n",celltype,":",dtxsid,":",sid,"\n",suse,":",star)
        boxplot(y~xnew,main=main,cex.main=1,
                ylim=c(0.01,100),log="x",xlab="BMD (uM)",ylab="",
                horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),names=newnames,col=cols)
        for(bmd in c(100,10,1,0.1,0.001)) lines(c(bmd,bmd),c(0,1000000))
        lines(c(burst_pod,burst_pod),c(0,1000000),lwd=4,col="cyan")
        row[1,"name"] = chem
        row[1,"celltype"] = celltype
        row[1,"dtxsid"] = dtxsid
        row[1,"sid"] = sid
        row[1,"use"] = suse
        row[1,"annotation"] = star
        row[1,"pod"] = pod
        row[1,"burst_pod"] = burst_pod
        row[1,"ratio"] = burst_pod/pod
        row[1,"specific_targets"] = targets
        row[1,"pod_targets"] = pod_target
        summary = rbind(summary,row)
        if(!to.file) browser()
      }
    }
  }
  if(to.file) dev.off()
  file <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  write.xlsx(summary,file)
}

