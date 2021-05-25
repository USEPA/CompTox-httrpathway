#--------------------------------------------------------------------------------------
#' Generate chemical-wise boxplot of the BMD distributions by super_target
#'
#' @param to.file If TRUE, send the plots to a file
#' @param do.load If TRUE, load hte large HTTr data set into memory
#' @param dataset Name of the HTTr data set
#' @param sigcatalog Name of the signature catalog to use
#' @param sigset Name of the signature set
#' @param method Scoring method
#' @param celltype Name of the cell type
#' @param hccut Exclude rows in the data set with hitcall less than this value
#' @param tccut Exclude rows in the data set with top_over_cutoff less than this value
#' @param cutoff The minimum number of signatures hat have to be active in a super
#' target for the super target to be considered active
#' @param minconc Minimum concentration used in the plots
#' @param maxconc Maximum concentration used in the plots
#'
#' After running this, run the following ...
#' superTargetPODplot
#' superTargetStats
#--------------------------------------------------------------------------------------
superTargetBoxplot <- function(to.file=T,
                               do.load=T,
                               dataset="u2os_toxcast_pfas_pe1_normal_refchems",
                               sigcatalog="signatureDB_master_catalog 2021-04-24",
                               sigset="screen_large",
                               method="fc",
                               celltype="U2OS",
                               hccut=0.95,
                               tccut=1.5,
                               cutoff=5,
                               minconc=0.001,
                               maxconc=100) {
  printCurrentFunction(paste(dataset,sigset,method))

  if(!exists("TOXCAST")) {
    file = "../input/toxcast_matrix/toxcast_active_by_source.RData"
    load(file=file)
    TOXCAST <<- TOXCAST
  }
  toxcast = TOXCAST

  catalog = read.xlsx(paste0("../input/signatures/",sigcatalog,".xlsx"))
  catalog = catalog[catalog[,sigset]==1,]
  temp = unique(catalog[,c("parent","super_target")])
  sttot = table(temp$super_target)
  file = "../input/chemicals/PFAS synonyms.xlsx"
  synonyms = read.xlsx(file)
  synonyms = synonyms[!is.na(synonyms$nickname),]
  rownames(synonyms) = synonyms$dtxsid

  source.list = c("ACEA","ACEA_Cytotoxicity","APR_Cytotoxicity","APR_dn","APR_up",
                  "ATG_CIS","ATG_TRANS","BSK_Cytotoxicity","BSK_down","BSK_up","ZF_terata")

  toxcast.cytotox = TOXCAST[is.element(TOXCAST$source,c("ACEA_Cytotoxicity","APR_Cytotoxicity","BSK_Cytotoxicity")),]
  toxcast.bsk = TOXCAST[is.element(TOXCAST$source,c("BSK_down","BSK_up")),]
  toxcast.atg = TOXCAST[is.element(TOXCAST$source,c("ATG_CIS","ATG_TRANS")),]
  toxcast.other = TOXCAST[is.element(TOXCAST$source,c("ACEA","APR_dn","APR_up","ZF_terata")),]

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }
  mat = MAT

  dlist1 = unique(toxcast$dtxsid)
  dlist2 = unique(mat$dtxsid)
  cat("ToxCast chemicals:",length(dlist1),"\n")
  cat("HTTr chemicals:",length(dlist2),"\n")

  mat[is.na(mat$bmd),"bmd"] = 1000
  mat[is.na(mat$top_over_cutoff),"top_over_cutoff"] = 0
  mat[mat$hitcall<hccut,"bmd"] = 1000
  mat[mat$hitcall<hccut,"hitcall"] = 0
  mat[mat$top_over_cutoff<tccut,"bmd"] = 1000
  mat[mat$top_over_cutoff<tccut,"hitcall"] = 0

  st.list = sort(unique(mat$super_target))
  name.list = c("sample_id","dtxsid","name","use_class","targets","super_target","chem_super_target","bmd_median","match_chem","active","nhit_up","nhit_dn")
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
      #print(chem_super_target)
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
      res1chem$nhit_up = 0
      res1chem$nhit_dn = 0
      res1chem$newname = "-"
      for(st in st.list) {
        tempst = temp[is.element(temp$super_target,st),]
        tempst0 = tempst[tempst$hitcall==0,]
        tempst1 = tempst[tempst$hitcall>0,]
        if(nrow(tempst1)>=cutoff) {
          res1chem[is.element(res1chem$super_target,st),"nhit_up"] = nrow(tempst1[tempst1$top>0,])
          res1chem[is.element(res1chem$super_target,st),"nhit_dn"] = nrow(tempst1[tempst1$top<0,])
          res1chem[is.element(res1chem$super_target,st),"active"] = 1
          res1chem[is.element(res1chem$super_target,st),"bmd_median"] = median(tempst1$bmd)
        }
      }

      ######################################################################
      # deal with the hits
      ######################################################################
      #cat("hits\n")
      temp = temp1
      cat(chem,nrow(temp),"\n")
      if(nrow(temp)>0) {
        x = temp$super_target
        y = temp$bmd

        res0 = res1chem[res1chem$active==0,]
        res = res1chem[res1chem$active==1,]

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
        #cat(length(unique(xnew)),"\n")
        #browser()
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
          qb = quantile(bmds,probs=seq(0,1,0.05))
          pod_sig = qb[2]
          nsig = nrow(temp)
        }

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
        yy4 = NULL
        xx4 = NULL
        yy5 = NULL
        xx5 = NULL
        yy6 = NULL
        xx6 = NULL

        if(is.element(dtxsid,htpp$dtxsid)) {
          rtemp = htpp[is.element(htpp$dtxsid,dtxsid),]
          mcf7 = rtemp[is.element(rtemp$celltype,"MCF7"),"bmc"]
          mcf7 = mcf7[!is.na(mcf7)]
          u2os = rtemp[is.element(rtemp$celltype,"U2OS"),"bmc"]
          u2os = u2os[!is.na(u2os)]

          if(length(u2os)>0) {
            yy1 = u2os
            xx1 = u2os
            xx1[] = "990 HTPP U2OS"
          }
          if(length(mcf7)>0) {
            yy2 = mcf7
            xx2 = mcf7
            xx2[] = "991 HTPP MCF7"
          }
        }

        if(is.element(dtxsid,toxcast.cytotox$dtxsid)) {
          rtemp = toxcast[is.element(toxcast.cytotox$dtxsid,dtxsid),]
          tpod = rtemp$ac50
          if(length(tpod)>0) {
            yy3 = tpod
            xx3 = tpod
            xx3[] = "992 ToxCast Cytotox"
          }
        }
        if(is.element(dtxsid,toxcast.bsk$dtxsid)) {
          rtemp = toxcast[is.element(toxcast.bsk$dtxsid,dtxsid),]
          tpod = rtemp$ac50
          if(length(tpod)>0) {
            yy4 = tpod
            xx4 = tpod
            xx4[] = "993 ToxCast BSK"
          }
        }
        if(is.element(dtxsid,toxcast.atg$dtxsid)) {
          rtemp = toxcast[is.element(toxcast.atg$dtxsid,dtxsid),]
          tpod = rtemp$ac50
          if(length(tpod)>0) {
            yy5 = tpod
            xx5 = tpod
            xx5[] = "994 ToxCast ATG"
          }
        }
        if(is.element(dtxsid,toxcast.other$dtxsid)) {
          rtemp = toxcast[is.element(toxcast.other$dtxsid,dtxsid),]
          tpod = rtemp$ac50
          if(length(tpod)>0) {
            yy6 = tpod
            xx6 = tpod
            xx6[] = "995 ToxCast Other"
          }
        }
         if(length(y)>0 && nrow(res)>1) {
          newnames = res$super_target
          #cat("[1] ",length(unique(xnew)),length(newnames),"\n")
          #browser()

          cols = newnames
          cols[] = "white"
          for(k in 1:length(newnames)) {
            nnk = newnames[k]
            nnk0 = nnk
            top = temp[is.element(temp$super_target,nnk),"top"]
            if(length(top)>0) {
              np = length(top[top>0])
              nm = length(top[top<0])
              ntot = sttot[nnk]
              nnk = paste0(nnk," (",np,",",nm,"|",ntot,")")
            }
            if(is.element(nnk0,chem_super_target)) {
              nnk = paste(nnk,"*")
              cols[k] = "red"
            }
            if(is.element(nnk0,"Stress")) {
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
            newnames = c(newnames,"ToxCast Cytotox")
            cols = c(cols,"orange")
          }
          if(length(yy4)>0) {
            y = c(y,yy4)
            xnew = c(xnew,xx4)
            newnames = c(newnames,"ToxCast BSK")
            cols = c(cols,"white")
          }
          if(length(yy5)>0) {
            y = c(y,yy5)
            xnew = c(xnew,xx5)
            newnames = c(newnames,"ToxCast ATG")
            cols = c(cols,"white")
          }
          if(length(yy6)>0) {
            y = c(y,yy6)
            xnew = c(xnew,xx6)
            newnames = c(newnames,"ToxCast Other")
            cols = c(cols,"white")
          }
          #cat("[2] ",length(unique(xnew)),length(newnames),"\n")

          if(is.element(dtxsid,synonyms$dtxsid)) {
            nn = synonyms[dtxsid,"nickname"]
            chem = paste(nn,":",chem)
          }
          main=paste(chem,"\n",celltype,":",dtxsid,"\n",sid,"\n",suse,"\n",star)
          #browser()
          boxplot(y~xnew,main=main,cex.main=0.8,
                  ylim=c(minconc,100),log="x",xlab="BMD (uM)",ylab="",
                  horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),names=newnames,col=cols)
          for(bmd in c(100,10,1,0.1,0.01,0.001,0.0001,0.00001)) lines(c(bmd,bmd),c(0,1000000))
          median.signature.bmd
          lines(c(burst_pod_sig,burst_pod_sig),c(0,1000000),lwd=4,col="cyan")
          #lines(c(burst_pod_st,burst_pod_st),c(0,1000000),lwd=4,col="cyan")
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

