#--------------------------------------------------------------------------------------
#' Generate chemicalwise boxplot of the BMD distributions by super_target
#'
#'  heparg2d_toxcast_pfas_pe1_normal
#'  mcf7_ph1_pe1_normal_block_123
#'  u2os_toxcast_pfas_pe1_normal
#'  PFAS_HepaRG
#'  PFAS_U2OS
#'  u2os_pilot_pe1_normal_null_pilot_lowconc
#'  u2os_toxcast_pfas_pe1_normal_refchems
#'  heparg2d_toxcast_pfas_pe1_normal_refchems
#'
#' After running this, run the following ...
#' superTargetPODplot
#' superTargetStats
#--------------------------------------------------------------------------------------
superTargetBoxplot <- function(to.file=F,
                               do.load=F,
                               dataset="PFAS_U2OS",
                               sigset="screen_large",
                               method="fc",
                               celltype="U2OS",
                               hccut=0.95,
                               tccut=1.5,
                               minconc=0.001,
                               maxconc=100) {
  printCurrentFunction(paste(dataset,sigset,method))

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
  #x1 = mat[is.element(mat$super_target,"RAR"),]
  #x1 = x1[is.element(x1$name,"all-trans-Retinoic acid"),]
  mat[is.na(mat$bmd),"bmd"] = 1000
  mat[is.na(mat$top_over_cutoff),"top_over_cutoff"] = 0
  mat[mat$hitcall<hccut,"bmd"] = 1000
  mat[mat$hitcall<hccut,"hitcall"] = 0
  mat[mat$top_over_cutoff<tccut,"bmd"] = 1000
  mat[mat$top_over_cutoff<tccut,"hitcall"] = 0
  #x2 = mat[is.element(mat$super_target,"RAR"),]
  #x2 = x2[is.element(x2$name,"all-trans-Retinoic acid"),]
  #browser()
  st.list = sort(unique(mat$super_target))
  name.list = c("sample_id","dtxsid","name","use_class","targets","super_target","chem_super_target","bmd_median","match_chem","active")
  res1chem = as.data.frame(matrix(nrow=length(st.list),ncol=length(name.list)))
  names(res1chem) = name.list

  file = "../htpp/htpp_category.xlsx"
  htpp = read.xlsx(file)

  x = mat$bmd
  x[x<minconc] = minconc
  mat$bmd = x

  file = "../input/chemicals/refchem_super_target_map.xlsx"
  stdb = read.xlsx(file)
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
  res.all = NULL
  #mat = mat[is.element(mat$name,"Butylparaben"),]
  chem.list = sort(unique(mat$name))
  for(chem in chem.list) {
    tempchem = mat[is.element(mat$name,chem),]
    dtxsid = tempchem[1,"dtxsid"]
    name = tempchem[1,"name"]
    chem_annotation = NULL
    chem_use = NULL
    chem_super_target = NULL

    if(is.element(dtxsid,stdb$dtxsid)) {
      info = stdb[is.element(stdb$dtxsid,dtxsid),]
      info = info[1,]
      chem_annotation = info[1,"target_all_information"]
      if(!is.na(info[1,"target_gene_family"])) chem_annotation = paste0(chem_annotation,":",info[1,"target_gene_family"])
      if(!is.na(info[1,"target_genes"])) chem_annotation = paste0(chem_annotation,":",info[1,"target_genes"])
      chem_use = info[1,"use_class"]
      chem_super_target = info[1,"stlist"]
      chem_super_target = str_split(chem_super_target,"\\|")[[1]]
      print(chem_super_target)
    }
    sid.list = unique(tempchem$sample_id)
    for(sid in sid.list) {
      temp = tempchem[is.element(tempchem$sample_id,sid),]
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

      temp1 = temp[temp$hitcall>0,]
      temp0 = temp[temp$hitcall==0,]
      #name.list = c("sample_id","dtxsid","name","use_class","targets","super_target",
      #              "chem_super_target","bmd_median","match_chem","hit")

      res1chem[] = NA
      res1chem$sample_id = sid
      res1chem$dtxsid = dtxsid
      res1chem$name = name
      res1chem$use_class = suse
      res1chem$targets = star
      res1chem$super_target = st.list
      res1chem$chem_super_target = paste(chem_super_target,collapse="|")
      res1chem$match_chem = 0
      res1chem[is.element(res1chem$super_target,chem_super_target),"match_chem"] = 1
      res1chem$bmd_median = 1000
      res1chem$active = 0
      res1chem$count = 0
      res1chem$newname = "-"
      for(st in st.list) {
        tempst = temp[is.element(temp$super_target,st),]
        tempst0 = tempst[tempst$hitcall==0,]
        tempst1 = tempst[tempst$hitcall>0,]
        if(nrow(tempst1)>0) {
          res1chem[is.element(res1chem$super_target,st),"active"] = 1
          res1chem[is.element(res1chem$super_target,st),"bmd_median"] = median(tempst1$bmd)
        }
      }
      #if(nrow(temp1)>0) {
      #  x = temp1$super_target
      #  y = temp1$bmd
      #  for(i in 1:length(st.list)) {
      #    st = st.list[i]
      #    rtemp = y[is.element(x,st)]
      #    res[i,"bmd_median"] = median(rtemp)
      #  }
      #}

      ######################################################################
      # deal with the hits
      ######################################################################
      cat("hits\n")
      temp = temp1
      cat(chem,nrow(temp),"\n")
      if(nrow(temp)>0) {
        x = temp$super_target
        y = temp$bmd

        #st.list = unique(x)
        #name.list = c("sample_id","dtxsid","name","use_class","targets","super_target","stlist","bmd_median","match_chem")
        #res = as.data.frame(matrix(nrow=length(st.list),ncol=length(name.list)))
        #names(res) = name.list
        #res$sample_id = sid
        #res$dtxsid = dtxsid
        #res$name = name
        #res$super_target = st.list
        #res$match_chem = 0
        #res$use_class = suse
        #res$targets = star
        #res$stlist = chem_super_target
        #res[is.element(res$super_target,chem_super_target),"match_chem"] = 1
        #for(i in 1:length(st.list)) {
        #  st = st.list[i]
        #  rtemp = y[is.element(x,st)]
        #  res[i,"bmd_median"] = median(rtemp)
        #}
        #res$active = 1
        #res$newname = ""
        res0 = res1chem[res1chem$active==0,]
        res = res1chem[res1chem$active==1,]

        #cat("1:",length(st.list),nrow(res0),nrow(res),"\n")
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
        #cat("2:",length(st.list),nrow(res0),nrow(res),"\n")

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
        res.all = rbind(res.all,res,res0)
        #cat("3:",length(st.list),nrow(res0),nrow(res),"\n")

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
        #cat("4:",length(st.list),nrow(res0),nrow(res),"\n")

        if(nrow(temp)>2) {
          bmds = sort(temp$bmd)
          db = density(bmds)
          burst_pod_sig = db$x[which.max(db$y)]
          #burst_pod_sig = 10**(median(log10(temp0$bmd)))
          qb = quantile(bmds,probs=seq(0,1,0.05))
          pod_sig = qb[2]
          nsig = nrow(temp)
        }
        #cat("5:",length(st.list),nrow(res0),nrow(res),"\n")

         nmax = 40
        if(nrow(res)>nmax) {
          res$useme = 0
          res[res$match_chem==1,"useme"] = 1
          res[is.element(res$super_target,"Stress"),"useme"] = 1
          n0 = nrow(res[res$useme==1,])
          for(kk in 1:nrow(res)) {
            if(n0<nmax) {
              if(res[kk,"useme"]==0) {
                res[kk,"useme"]=1
                n0 = n0 +1
              }
            }
          }

          mask = xnew
          mask[] = 0
          res = res[res$useme==1,]
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
            if(is.element(nnk,chem_super_target)) {
              nnk = paste(nnk,"*")
              cols[k] = "red"
            }
            if(is.element(nnk,"Stress")) {
              cols[k] = "orange"
            }
            newnames[k] = nnk
          }
          if(length(yy1)>0) {
            y = c(y,yy1)
            xnew = c(xnew,xx1)
            newnames = c(newnames,"HTPP U2OS")
            cols = c(cols,"white")
          }
          if(length(yy2)>0) {
            y = c(y,yy2)
            xnew = c(xnew,xx2)
            newnames = c(newnames,"HTPP MCF7")
            cols = c(cols,"white")
          }
          if(length(yy3)>0) {
            y = c(y,yy3)
            xnew = c(xnew,xx3)
            newnames = c(newnames,"ToxCast")
            cols = c(cols,"white")
          }
          main=paste(chem,"\n",celltype,":",dtxsid,"\n",sid,"\n",suse,"\n",star)
          boxplot(y~xnew,main=main,cex.main=0.8,
                  ylim=c(minconc,100),log="x",xlab="BMD (uM)",ylab="",
                  horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),names=newnames,col=cols)
          for(bmd in c(100,10,1,0.1,0.01,0.001,0.0001,0.00001)) lines(c(bmd,bmd),c(0,1000000))
          lines(c(burst_pod_st,burst_pod_st),c(0,1000000),lwd=4,col="cyan")
          counter = counter+1
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
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_summary.xlsx")
  write.xlsx(summary,file)
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_all.RData")
  save(res.all,file=file)
  cat(counter,nrow(summary),"\n")
}

