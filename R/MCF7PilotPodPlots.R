#' Carry out analyses of different POD methods for the pilot study
#'
MCF7PilotPodPlots <- function(to.file=F,
                              method="gsea",
                              celltype="MCF7",
                              sigset="screen_large",
                              hccut=0.9,
                              tccut=1,
                              cutoff=3,
                              condition="all") {
  printCurrentFunction()

  ######################################################################################################
  # read the httrpathway data
  ######################################################################################################
  nset = 6
  dataset.list = c(
    "MCF7_pilot_DMEM_6hr_pilot_normal_pe_1",
    "MCF7_pilot_DMEM_12hr_pilot_normal_pe_1",
    "MCF7_pilot_DMEM_24hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_6hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_12hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_24hr_pilot_normal_pe_1"
  )
  media.list = c("DMEM","DMEM","DMEM","PRF","PRF","PRF")
  time.list = c(6,12,24,6,12,24)
  res = NULL
  for(i in 1:nset) {
    dataset = dataset.list[i]
    file = paste0("../output/signature_pod/signature_pod_",condition,"_",sigset,"_",dataset,"_",method,"_",hccut,"_",cutoff,".xlsx")
    temp = read.xlsx(file)
    temp$media = media.list[i]
    temp$time = time.list[i]
    res = rbind(res,temp)
  }
  res$condition = paste(res$media,res$time)
  nlist = names(res)
  nlist[is.element(nlist,"signature_pod_95")] = "sig_05"
  nlist[is.element(nlist,"signature_pod_95.lci")] = "sig_05.lci"
  nlist[is.element(nlist,"signature_pod_95.uci")] = "sig_05.uci"
  nlist[is.element(nlist,"gene_pod_95")] = "gene_05"
  nlist[is.element(nlist,"gene_pod_95.lci")] = "gene_05.lci"
  nlist[is.element(nlist,"gene_pod_95.uci")] = "gene_05.uci"
  nlist[is.element(nlist,"super_target_pod")] = "sig_target"

  nlist[is.element(nlist,"signature_pod_min")] = "sig_min"
  nlist[is.element(nlist,"signature_pod_min.lci")] = "sig_min.lci"
  nlist[is.element(nlist,"signature_pod_min.uci")] = "sig_min.uci"
  nlist[is.element(nlist,"signature_pod_abs5")] = "sig_abs5"
  nlist[is.element(nlist,"signature_pod_abs5.lci")] = "sig_abs5.lci"
  nlist[is.element(nlist,"signature_pod_abs5.uci")] = "sig_abs5.uci"
  nlist[is.element(nlist,"gene_pod_min")] = "gene_min"
  nlist[is.element(nlist,"gene_pod_min.lci")] = "gene_min.lci"
  nlist[is.element(nlist,"gene_pod_min.uci")] = "sig_min.uci"
  nlist[is.element(nlist,"gene_pod_abs5")] = "gene_abs5"
  nlist[is.element(nlist,"gene_pod_abs5.lci")] = "gene_abs5.lci"
  nlist[is.element(nlist,"gene_pod_abs5.uci")] = "gene_abs5.uci"
  nlist[is.element(nlist,"gene_burst_pod")] = "gene_burst"
  names(res) = nlist

  ######################################################################################################
  # read the httrpathway / bmdexpress form data
  ######################################################################################################
  dir = "../output/mcf7_pilot/"

  # file = paste0(dir,"tcplfit2_per_chemical_BMDs_9_8_2021.rds")
  if(!exists("BMDS2")) {
    cat("read bmds2 data file\n")
    file = paste0(dir,"tcplfit2_BMDs_9_7_2021.rds")
    bmds2 = readRDS(file=file)
    BMDS2 <<- bmds2
  }
  bmds2 = BMDS2
  res$bmds2_gene_min = 1000
  res$bmds2_gene_abs5 = 1000
  res$bmds2_gene_05 = 1000
  for(i in 1:nrow(res)) {
    dtxsid = res[i,"dtxsid"]
    media = res[i,"media"]
    time = res[i,"time"]
    temp1 = bmds2[is.element(bmds2$chemical,dtxsid),]
    condition = paste0("GENE_CR_MCF7_pilot_",media,"_",time,"hr_pilot_normal_pe_1_0")
    temp2 = temp1[is.element(temp1$condition,condition),]
    bmds = sort(as.numeric(temp2$bmd))
    q = quantile(bmds,probs=seq(0,1,0.05))
    res[i,"bmds2_gene_min"] = bmds[1]
    res[i,"bmds2_gene_abs5"] = bmds[5]
    res[i,"bmds2_gene_05"] = q[2]
  }
  ######################################################################################################
  # read the bmdExpress data
  ######################################################################################################
  file = paste0(dir,"bmdExpress_mcf7_pilot_full_summary.xlsx")
  bmds = read.xlsx(file)
  bmds[is.element(bmds$media,"PRFDMEM"),"media"] = "PRF"
  res$bmds_gene_min = 1000
  res$bmds_gene_05 = 1000
  res$bmds_sig_min = 1000
  res$bmds_sig_05 = 1000
  for(i in 1:nrow(res)) {
    dtxsid = res[i,"dtxsid"]
    media = res[i,"media"]
    time = res[i,"time"]
    temp1 = bmds[is.element(bmds$dtxsid,dtxsid),]
    temp2 = temp1[is.element(temp1$timeh,time),]
    temp3 = temp2[is.element(temp2$media,media),]
    res[i,"bmds_gene_min"] = temp3[1,"min_gene_bmd"]
    res[i,"bmds_gene_05"] = temp3[1,"q05_gene_bmd"]
    res[i,"bmds_sig_min"] = temp3[1,"min_sig_bmd"]
    res[i,"bmds_sig_05"] = temp3[1,"q05_sig_bmd"]
  }

  btypes = c("bmds_gene_min","bmds_gene_05","bmds_sig_min","bmds_sig_05")
  for(btype in btypes) {
    x = res[,btype]
    x[x<0.003] = 0.003
    res[,btype] = x
  }
  ######################################################################################################
  # read the plier data
  ######################################################################################################
  file = paste0(dir,"plier_lv_pod.xlsx")
  plier = read.xlsx(file)
  res$plier_lv_min = NA
  res$plier_lv_abs5 = NA
  res$plier_lv_05 = NA
  for(i in 1:nrow(res)) {
    dtxsid = res[i,"dtxsid"]
    condition = res[i,"condition"]
    temp = plier[is.element(plier$dtxsid,dtxsid),]
    temp = temp[is.element(temp$condition,condition),]
    res[i,"plier_lv_min"] = temp[1,"plier_lv_pod_min"]
    res[i,"plier_lv_abs5"] = temp[1,"plier_lv_pod_abs5"]
    res[i,"plier_lv_05"] = temp[1,"plier_lv_pod_95"]
  }

   if(to.file) {
    fname = paste0(dir,"mcf7_pilot_pod_",method,"_",hccut,"_",tccut,"_",cutoff,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,5,6,2))

  tlist = c("gene_min",       "gene_05",          "gene_abs5",
            "sig_min",        "sig_05",           "sig_abs5",
            "sig_target",
            "plier_lv_min",   "plier_lv_05",      "plier_lv_abs5",
            "bmds_gene_min",  "bmds_gene_05",
            "bmds_sig_min",   "bmds_sig_05",
            "bmds2_gene_min", "bmds2_gene_abs5",  "bmds2_gene_05"
  )

  for(type in tlist) {
    x = res[,type]
    x[is.na(x)] = 1000
    x[x<0.003] = 0.003
    res[,type] = x
  }

  col.list = c("blue","blue","blue",
               "green","green","green",
               "red",
               "gray","gray","gray",
               "orange","orange",
               "black","black",
               "cyan","cyan","cyan")
  pch.list = c(25,21,22,
               25,21,22,
               25,
               25,21,22,
               25,21,
               25,21,
               25,22,21)
  cex.list = c(1,1,1,
               1,1,1,
               1,
               1,1,1,
               1,1,
               1,1,
               1,1,1)

  ####################################################################################
  # compare the PODs to the ToxCastvalues
  ####################################################################################
  file = paste0(dir,"toxcast_pods_wide_invitrodb_v3_4_09072021.xlsx")
  toxcast = read.xlsx(file)
  rownames(toxcast) = toxcast$dsstox_substance_id
  dlist = res$dtxsid
  nchem = length(dlist)
  method = NULL
  cdiff = NULL
  for(k in 1:6) {
    media = media.list[k]
    time = time.list[k]
    temp = res[is.element(res$media,media),]
    temp = temp[temp$time==time,]
    rownames(temp) = temp$dtxsid
    # plot(c(1,1),type="n",log="x",xlim=c(7e-6,1e2),ylim=c(0,nchem),xlab="pod (uM)",
    #      ylab="",main=paste("ER Actives\n",media," ",time,"hr"),cex.lab=1.4,cex.axis=1.4,
    #      yaxt="n",xaxt="n")
    # axis(side=1,at=c(1e-3,1e-2,1e-1,1,10,100),cex.lab=1.4,cex.axis=1.4)
    for(i in 1:nchem) {
      y = i-0.5
      dtxsid = dlist[i]
      # name = ermodel[dtxsid,"name"]
      # text(1e-3,y,name,pos=2,cex=1.5)
      # lines(c(1e-10,1e10),c(i,i))
      pod0 = toxcast[dtxsid,"toxcast_min_modl_acc_top12"]
      pod0 = 10**pod0
      # lines(c(pod0,pod0),c(i-1,i),col="red",lwd=2)
      for(j in 1:length(tlist)) {
        type = tlist[j]
        pod = temp[dtxsid,type]
        # points(pod,y,pch=pch.list[j],bg=col.list[j],cex=cex.list[j]*2)
        method = c(method,type)
        delta = log10(pod) - log10(pod0)
        cdiff = c(cdiff,delta)
      }
    }
  }

  par(mfrow=c(1,1),mar=c(4,15,3,2))
  boxplot(cdiff~method,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),
          main="Prediction Ratio for ToxCast",xlab="log POD ratio",ylab="",ylim=c(-3,8),cex.axis=1.3,cex.lab=1.3)
  for(x in seq(from=-3,to=8)) lines(c(x,x),c(0,100),col="gray")
  lines(c(0,0),c(0,100),col="red",lwd=2)
  lines(c(-1,-1),c(0,100),col="red",lwd=2)
  lines(c(1,1),c(0,100),col="red",lwd=2)
  if(!to.file) browser()

  ####################################################################################
  # compare the PODs to the ER model values
  ####################################################################################
  file = paste0(dir,"ER_chems all mcf7_ph1_pe1_normal_block_123_allPG screen_large 0.9 10.xlsx")
  ermodel = read.xlsx(file)
  erchems = c(
    "Fulvestrant",
    "4-Hydroxytamoxifen",
    "Clomiphene citrate (1:1)",
    "Bisphenol B",
	  "Bisphenol A",
    "4-Nonylphenol, branched",
    "4-Cumylphenol"
  )
  ermodel = ermodel[is.element(ermodel$name,erchems),]
  rownames(ermodel) = ermodel$name
  ermodel = ermodel[erchems,]
  rownames(ermodel) = ermodel$dtxsid
  nchem = nrow(ermodel)
  dlist = ermodel$dtxsid
  par(mfrow=c(3,1),mar=c(4,5,5,2))

  method = NULL
  cdiff = NULL
  for(k in 1:6) {
    media = media.list[k]
    time = time.list[k]
    temp = res[is.element(res$media,media),]
    temp = temp[temp$time==time,]
    rownames(temp) = temp$dtxsid
    plot(c(1,1),type="n",log="x",xlim=c(7e-6,1e2),ylim=c(0,nchem),xlab="pod (uM)",
         ylab="",main=paste("ER Actives\n",media," ",time,"hr"),cex.lab=1.4,cex.axis=1.4,
         yaxt="n",xaxt="n")
    axis(side=1,at=c(1e-3,1e-2,1e-1,1,10,100),cex.lab=1.4,cex.axis=1.4)
    for(i in 1:nchem) {
      y = i-0.5
      dtxsid = dlist[i]
      name = ermodel[dtxsid,"name"]
      text(1e-3,y,name,pos=2,cex=1.5)
      lines(c(1e-10,1e10),c(i,i))
      pod0 = min(ermodel[dtxsid,"hts.pod.agonist"],ermodel[dtxsid,"hts.pod.antagonist"])
      pod0 = 10**pod0
      lines(c(pod0,pod0),c(i-1,i),col="red",lwd=2)
      for(j in 1:length(tlist)) {
        type = tlist[j]
        pod = temp[dtxsid,type]
        points(pod,y,pch=pch.list[j],bg=col.list[j],cex=cex.list[j]*2)
        method = c(method,type)
        delta = log10(pod) - log10(pod0)
        cdiff = c(cdiff,delta)
      }
    }
  }
  if(!to.file) browser()
  par(mfrow=c(1,1),mar=c(4,15,3,2))
  boxplot(cdiff~method,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),
          main="Prediction Ratio for ER",xlab="log POD ratio",ylab="",ylim=c(-3,3),cex.axis=1.3,cex.lab=1.3)
  lines(c(0,0),c(0,100),col="red",lwd=2)
  lines(c(-1,-1),c(0,100),col="red",lwd=2)
  lines(c(1,1),c(0,100),col="red",lwd=2)
  if(!to.file) browser()

  ####################################################################################
  # compare the PODs to the statin values
  ####################################################################################
  schems = c(
    "Simvastatin",
    "Lovastatin"
  )
  spods = c(11.2,8.4)
  # references
  # simvastatin https://onlinelibrary.wiley.com/doi/pdf/10.1002/clc.4960261507
  # lovastatin https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1742-7843.2005.pto_134.x
  nchem = nrow(schems)
  statins = unique(res[is.element(res$name,schems),c("dtxsid","name")])
  rownames(statins) = statins$dtxsid
  statins$pod0 = NA
  statins[is.element(statins$name,"Simvastatin"),"pod0"] = 11.2
  statins[is.element(statins$name,"Lovastatin"),"pod0"] = 8.4
  dlist = statins$dtxsid
  nchem = length(dlist)
  par(mfrow=c(3,1),mar=c(4,5,5,2))

  method = NULL
  cdiff = NULL
  for(k in 1:6) {
    media = media.list[k]
    time = time.list[k]
    temp = res[is.element(res$media,media),]
    temp = temp[temp$time==time,]
    rownames(temp) = temp$dtxsid
    plot(c(1,1),type="n",log="x",xlim=c(1e-3,1e2),ylim=c(0,nchem),xlab="pod (uM)",
         yaxt="n",xaxt="n",
         ylab="",main=paste("HMGCR Actives\n",media," ",time,"hr"),cex.lab=1.2,cex.axis=1.2)
    axis(side=1,at=c(0.01,0.1,1,10,100))
    for(i in 1:nchem) {
      y = i-0.5
      dtxsid = dlist[i]
      name = statins[dtxsid,"name"]
      text(1e-2,y,name,pos=2,cex=1.5)
      lines(c(1e-10,1e10),c(i,i))
      pod0 = statins[dtxsid,"pod0"]
      lines(c(pod0,pod0),c(i-1,i),col="red",lwd=2)
      for(j in 1:length(tlist)) {
        type = tlist[j]
        pod = temp[dtxsid,type]
        points(pod,y,pch=pch.list[j],bg=col.list[j],cex=cex.list[j]*2)
        method = c(method,type)
        delta = log10(pod) - log10(pod0)
        cdiff = c(cdiff,delta)
      }
    }
  }
  if(!to.file) browser()
  par(mfrow=c(1,1),mar=c(4,15,3,2))
  boxplot(cdiff~method,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),
          main="Prediction Ratio for HMGCR",xlab="log POD ratio",ylab="",ylim=c(-3,3))
  lines(c(0,0),c(0,100),col="red",lwd=2)
  lines(c(-1,-1),c(0,100),col="red",lwd=2)
  lines(c(1,1),c(0,100),col="red",lwd=2)
  if(!to.file) browser()

  ####################################################################################
  # deltas - fix everything and look at metrics
  ####################################################################################
  par(mfrow=c(5,4),mar=c(4,6,3,1))
  type0 = "sig_05"
  for(type in tlist) {
    if(type!=type0) {
      x = NULL
      y = NULL
      for(time in c(6,12,24)) {
        temp1 = res[res$time==time,]
        for(media in c("DMEM","PRF")) {
          temp2 = temp1[is.element(temp1$media,media),]
          for( time in c(6,12,24)) {
            temp3 = temp2[temp2$time==time,]
            for(dtxsid in unique(res$dtxsid)) {
              temp4 = temp3[is.element(temp3$dtxsid,dtxsid),]
              value = log10(temp4[1,type]/temp4[1,type0])
              tag = paste0(media," t=",time)
              x = c(x,tag)
              y = c(y,value)
            }
          }
        }
      }
      boxplot(y~x,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),main=type,
              xlab="log POD ratio",ylab="",ylim=c(-3,3))
      lines(c(0,0),c(0,100),col="red",lwd=2)
      lines(c(-1,-1),c(0,100),col="red",lwd=1)
      lines(c(1,1),c(0,100),col="red",lwd=1)
    }
  }
  if(!to.file) browser()
  ####################################################################################
  # deltas - fix everything and look at media
  ####################################################################################
  par(mfrow=c(5,4),mar=c(4,6,3,1))
  for(type in tlist) {
    x = NULL
    y = NULL
    for(time in c(6,12,24)) {
      temp1 = res[res$time==time,]
      for(dtxsid in unique(res$dtxsid)) {
        temp2 = temp1[is.element(temp1$dtxsid,dtxsid),]
        value = log10(temp2[temp2$media=="DMEM",type]/temp2[temp2$media=="PRF",type])
        tag = paste0("DMEM/PRF t=",time)
        x = c(x,tag)
        y = c(y,value)
      }
    }
    boxplot(y~x,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),main=type,xlab="log POD ratio",
            ylab="",ylim=c(-3,3),cex.axis=0.9,cex.lab=0.9)
    lines(c(0,0),c(0,100),col="red",lwd=2)
  }
  if(!to.file) browser()
  ####################################################################################
  # deltas - fix everything and look at time
  ####################################################################################
  par(mfrow=c(5,4),mar=c(4,6,3,1))
  for(type in tlist) {
    x = NULL
    y = NULL
    for(media in unique(res$media)) {
      temp1 = res[is.element(res$media,media),]
      for(dtxsid in unique(res$dtxsid)) {
        temp2 = temp1[is.element(temp1$dtxsid,dtxsid),]
        value = log10(temp2[temp2$time==12,type]/temp2[temp2$time==6,type])
        tag = paste(media,"12/6")
        x = c(x,tag)
        y = c(y,value)
        value = log10(temp2[temp2$time==24,type]/temp2[temp2$time==6,type])
        tag = paste(media,"24/6")
        x = c(x,tag)
        y = c(y,value)
      }
    }
    boxplot(y~x,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),main=type,xlab="log POD ratio",ylab="",ylim=c(-3,3))
    lines(c(0,0),c(0,100),col="red",lwd=2)
  }
  if(!to.file) browser()

  ####################################################################################
  # by condition
  ####################################################################################
  par(mfrow=c(1,1),mar=c(4,5,6,2))
  for(condition in unique(res$condition)) {
    temp = res[is.element(res$condition,condition),]
    temp = temp[order(temp$sig_05),]
    ptemp = plier[is.element(plier$condition,condition),]
    rownames(ptemp) = ptemp$dtxsid

    plot(c(1,1),type="n",main=paste("By Condition\n",condition),xlim=c(1e-6,1000),ylim=c(1,44),
         cex.axis=1.2,cex.lab=1.2,xlab="POD (uM)",ylab="",log="x",yaxt="n",xaxt="n")
    axis(side=1,at=c(0.001,0.01,0.1,1,10,100,1000))
    for(x in c(0.001,0.01,0.1,1,10,100,1000)) lines(c(x,x),c(0,100),col="gray")
    for(i in 1:nrow(temp)) {
      name = temp[i,"name"]
      dtxsid = temp[i,"dtxsid"]
      lines(c(1e-6,1000),c(i,i),col="gray")
      text(1e-3,i+0.5,name,pos=2,cex=0.8)
      for(j in 1:length(tlist)) {
        type = tlist[j]
        #if(j<=7)
          pod = temp[i,type]
        #else pod = ptemp[dtxsid,type]
        points(pod,i+0.5,pch=pch.list[j],bg=col.list[j],cex=cex.list[j])
      }
    }
    if(!to.file) browser()
  }

  ###############################################################################
  par(mfrow=c(1,1),mar=c(4,5,6,2))
  clist = unique(res$condition)
  chems = unique(res[,c("dtxsid","name")])

  res2 = NULL
  for(type in tlist) {
    chems$type = type
    vals = as.data.frame(matrix(nrow=nrow(chems),ncol=length(clist)))
    names(vals) = clist
    for(condition in clist) {
      temp1 = res[is.element(res$condition,condition),]
      for(i in 1:nrow(chems)) {
        dtxsid = chems[i,"dtxsid"]
        pod = temp1[is.element(temp1$dtxsid,dtxsid),type]
        vals[i,condition] = pod
      }
    }
    x = cbind(chems,vals)
    res2 = rbind(res2,x)
  }
  temp = res2[,clist]
  res2$minval = apply(temp,FUN=min,MARGIN=1)

  col.list = c("yellow","orange","red","yellow","orange","red")
  pch.list = c(24,24,24,25,25,25)

  for(type in tlist) {
    temp = res2[is.element(res2$type,type),]

    temp = temp[order(temp$minval),]
    plot(c(1,1),type="n",main=paste("By Type\n",type),xlim=c(1e-6,1000),ylim=c(0,44),
         cex.axis=1.2,cex.lab=1.2,xlab="POD (uM)",ylab="",log="x",yaxt="n",xaxt="n")
    axis(side=1,at=c(0.001,0.01,0.1,1,10,100,1000))
    for(i in 1:nrow(temp)) {
      name = temp[i,"name"]
      dtxsid = temp[i,"dtxsid"]
      lines(c(1e-6,1000),c(i,i),col="gray")
      text(1e-3,i+0.5,name,pos=2,cex=0.8)
      #points(toxcast[dtxsid,"pod_uM"],i+0.5,pch=4)
      for(j in 1:length(clist)) {
        condition = clist[j]
        pod = temp[i,condition]
        points(pod,i+0.5,pch=pch.list[j],bg=col.list[j])
      }
    }
    if(!to.file) browser()
  }

  ###############################################################################
  par(mfrow=c(3,1),mar=c(4,20,4,8))
  x = NULL
  y = NULL
  for(condition in unique(res$condition)) {
    temp = res[is.element(res$condition,condition),]
    yy = NULL
    for(i in 1:nrow(temp)) {
      vals = log10(temp[i,tlist])
      maxval = max(vals,na.rm=T)
      minval = min(vals,na.rm=T)
      yy = c(yy,maxval-minval)
    }

    xx = yy
    condition = str_replace(condition," 6"," 06")
    xx[] = condition
    x = c(x,xx)
    y = c(y,yy)
  }
  y = 3-y
  boxplot(y~x,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),main="Range of PODs\nby Media and Time",
          xlab="log range",ylab="",ylim=c(-3,3))

  x = NULL
  y = NULL
  dlist = unique(res$dtxsid)
  for(type in tlist) {
    for(dtxsid in dlist) {
      temp = res[is.element(res$dtxsid,dtxsid),]
      logx = log10(temp[,type])
      maxval = max(logx)
      minval = min(logx)
      y = c(y,maxval-minval)
      x = c(x,type)
    }
  }
  y = 3-y
  boxplot(y~x,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),main="Range of PODs\nby POD Type",
          xlab="log range",ylab="",ylim=c(-3,3))
  for(x in c(-3,-2,-1,0,1,2,3)) lines(c(x,x),c(0,100),col="gray")

  if(!to.file) browser()
  if(to.file) dev.off()
}
