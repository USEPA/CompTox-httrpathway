#--------------------------------------------------------------------------------------
#' Generate chemicalwise boxplot of the BMD distributions by super_target
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
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
                               tccut=1.5) {
  printCurrentFunction(paste(dataset,sigset,method))

  #file = "../toxcast/hts.xlsx"
  #toxcast = read.xlsx(file)
  if(!exists("TOXCAST")) {
    load(file='../toxcast/toxcast_master.RData')
    temp = toxcast.master
    temp = temp[,c("dsstox_substance_id","modl_ga_uM")]
    names(temp) = c("dtxsid","pod")
    TOXCAST <<- temp
  }
  toxcast = TOXCAST
  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }
  mat = MAT
  mat[is.na(mat$bmd),"bmd"] = 1000
  mat[is.na(mat$top_over_cutoff),"top_over_cutoff"] = 0
  mat[mat$hitcall<hccut,"bmd"] = 1000
  mat[mat$hitcall<hccut,"hitcall"] = 0
  mat[mat$top_over_cutoff<tccut,"bmd"] = 1000
  mat[mat$top_over_cutoff<tccut,"hitcall"] = 0

  file = "../htpp/htpp_category.xlsx"
  htpp = read.xlsx(file)


  x = mat$bmd
  x[x<0.001] = 0.001
  mat$bmd = x
  file = "../input/signatures/refchemdb_chem_filtered_unique_with_new.xlsx"
  file = "../input/chemicals/httr_chemical_annotations 2020-09-14.xlsx"
  rcdb = read.xlsx(file)
  if(to.file) {
    fname <- paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,15,6,2))

  name.list = c("name","celltype","dtxsid","sid","use","annotation",
                "pod_st","burst_pod_st","nst",
                "pod_sig","burst_pod_sig","nsig",
                "ratio_st",
                "ratio_sig",
                "specific_targets",
                "pod_target")
  row = as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(row) = name.list
  summary = NULL
  counter = 0

  chem.list = sort(unique(mat$name))
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

      suse = paste(chem_use,collapse="|")
      star = paste(chem_annotation,collapse="|")
      pod_st = 1000
      burst_pod_st = 1000
      pod_sig = 1000
      burst_pod_sig = 1000
      nst = 0
      nsig = 0
      targets = "-"
      pod_target = "-"

      temp = temp[temp$hitcall>0,]

      cat(chem,nrow(temp),"\n")
      if(nrow(temp)>0) {
        x = temp$super_target
        y = temp$bmd

        st.list = unique(x)

        res = as.data.frame(matrix(nrow=length(st.list),ncol=2))
        names(res) = c("super_target","bmd_median")
        res$super_target = st.list
        for(i in 1:length(st.list)) {
          st = st.list[i]
          rtemp = y[is.element(x,st)]
          res[i,2] = median(rtemp)
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

        xnew = xnew[mask==1]
        y = y[mask==1]

        res = res[is.element(res$newname,xnew),]
        burst_pod_st = 1000
        pod_st = 1000
        burst_pod_sig = 1000
        pod_sig = 1000
        nst = 0
        nsig = 0
        targets = "-"
        pod_target = "-"
        if(nrow(res)>1) {
          burst_pod_st = 10**(median(log10(res$bmd_median)))
          pod_st = min(res$bmd_median)
          pod_target = res[1,"super_target"]
          rtemp = res[res$bmd_median<burst_pod_st/10,]
          targets = "-"
          if(nrow(rtemp)>0) targets = paste(rtemp$super_target,collapse="|")
          nst = nrow(res)
        }
        if(nrow(temp)>2) {
          bmds = sort(temp$bmd)
          db = density(bmds)
          burst_pod_sig = db$x[which.max(db$y)]
          #burst_pod_sig = 10**(median(log10(temp0$bmd)))
          qb = quantile(bmds,probs=seq(0,1,0.05))
          pod_sig = qb[2]
          nsig = nrow(temp)
          #browser()
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

        yy1 = NULL
        xx1 = NULL
        yy2 = NULL
        xx2 = NULL
        yy3 = NULL
        xx3 = NULL

        if(is.element(dtxsid,htpp$dtxsid)) {
          rtemp = htpp[is.element(htpp$dtxsid,dtxsid),]
          mcf7 = rtemp[is.element(rtemp$celltype,"MCF7"),"bmc"]
          mcf7 = mcf7[!is.na(mcf7)]
          u2os = rtemp[is.element(rtemp$celltype,"U2OS"),"bmc"]
          u2os = u2os[!is.na(u2os)]

          if(length(u2os)>0) {
            yy1 = u2os
            xx1 = u2os
            xx1[] = "997 HTPP U2OS"
          }
          if(length(mcf7)>0) {
            yy2 = mcf7
            xx2 = mcf7
            xx2[] = "998 HTPP MCF7"
          }
         }
        if(is.element(dtxsid,toxcast$dtxsid)) {
          rtemp = toxcast[is.element(toxcast$dtxsid,dtxsid),]
          #rtemp = rtemp[rtemp$hitcall==1,]
          tpod = rtemp$pod
          #tpod = tpod[!is.na(tpod)]
          #tpod = tpod[tpod<6]
          #tpod = 10**tpod

          if(length(tpod)>0) {
            yy3 = tpod
            xx3 = tpod
            xx3[] = "999 ToxCast"
          }
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
          if(length(yy1)>0) {
            y = c(y,yy1)
            xnew = c(xnew,xx1)
            newnames = c(newnames,"HTPP U2OS")
          }
          if(length(yy2)>0) {
            y = c(y,yy2)
            xnew = c(xnew,xx2)
            newnames = c(newnames,"HTPP MCF7")
          }
          if(length(yy3)>0) {
            y = c(y,yy3)
            xnew = c(xnew,xx3)
            newnames = c(newnames,"ToxCast")
          }
          main=paste(chem,"\n",celltype,":",dtxsid,":",sid,"\n",suse,":",star)
          boxplot(y~xnew,main=main,cex.main=1,
                  ylim=c(0.001,100),log="x",xlab="BMD (uM)",ylab="",
                  horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),names=newnames,col=cols)
          for(bmd in c(100,10,1,0.1,0.01,0.001)) lines(c(bmd,bmd),c(0,1000000))
          lines(c(burst_pod_st,burst_pod_st),c(0,1000000),lwd=4,col="cyan")
          counter = counter+1
          #if(is.element(dtxsid,toxcast$dtxsid)) {
          #  pod = toxcast[dtxsid,"pod_uM"]
          #  lines(c(pod,pod),c(0,1000000),lwd=4,col="red")
          #}

          if(!to.file) browser()
        }
      }

      row[1,"name"] = chem
      row[1,"celltype"] = celltype
      row[1,"dtxsid"] = dtxsid
      row[1,"sid"] = sid
      row[1,"use"] = suse
      row[1,"annotation"] = star
      row[1,"pod_st"] = pod_st
      row[1,"burst_pod_st"] = burst_pod_st
      row[1,"pod_sig"] = pod_sig
      row[1,"burst_pod_sig"] = burst_pod_sig
      row[1,"ratio_st"] = burst_pod_st/pod_st
      row[1,"ratio_sig"] = burst_pod_sig/pod_sig
      row[1,"nst"] = nst
      row[1,"nsig"] = nsig
      row[1,"specific_targets"] = targets
      row[1,"pod_target"] = pod_target
      summary = rbind(summary,row)

    }
  }
  if(to.file) dev.off()
  file <- paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".xlsx")
  write.xlsx(summary,file)
  cat(counter,nrow(summary),"\n")
}

