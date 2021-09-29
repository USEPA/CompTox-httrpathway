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
  dir = "../output/mcf7_pilot/"

  file = paste0(dir,"plier_lv_pod.xlsx")
  plier = read.xlsx(file)
  res$plier_lv_pod_min = NA
  res$plier_lv_pod_abs5 = NA
  res$plier_lv_pod_95 = NA
  for(i in 1:nrow(res)) {
    dtxsid = res[i,"dtxsid"]
    condition = res[i,"condition"]
    temp = plier[is.element(plier$dtxsid,dtxsid),]
    temp = temp[is.element(temp$condition,condition),]
    res[i,"plier_lv_pod_min"] = temp[1,"plier_lv_pod_min"]
    res[i,"plier_lv_pod_abs5"] = temp[1,"plier_lv_pod_abs5"]
    res[i,"plier_lv_pod_95"] = temp[1,"plier_lv_pod_95"]
  }

  if(to.file) {
    fname = paste0(dir,"mcf7_pilot_pod_",method,"_",hccut,"_",tccut,"_",cutoff,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,5,6,2))

  tlist = c("gene_pod_min","gene_pod_95","gene_pod_abs5",
            "signature_pod_min","signature_pod_95","signature_pod_abs5",
            "super_target_pod",
            "plier_lv_pod_min","plier_lv_pod_95","plier_lv_pod_abs5")
  col.list = c("blue","blue","blue","green","green","green","red","gray","gray","gray")
  pch.list = c(25,21,22,
               25,21,22,
               21,
               25,21,22)
  cex.list = c(0.0,1,0.0,
               0.5,1,0.5,
               1,
               0.5,1,0.5)

  ####################################################################################
  # deltas - fix everything and look at metrix
  ####################################################################################
  par(mfrow=c(5,2),mar=c(4,10,3,2))
  type0 = "signature_pod_95"
  for(type in tlist) {
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
    boxplot(y~x,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),main=type,xlab="log POD ratio",ylab="",ylim=c(-2,2))
    lines(c(0,0),c(0,100),col="red")
  }
  if(!to.file) browser()
  ####################################################################################
  # deltas - fix everything and look at media
  ####################################################################################
  par(mfrow=c(5,2),mar=c(4,10,3,2))
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
    boxplot(y~x,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),main=type,xlab="log POD ratio",ylab="",ylim=c(-2,2))
    lines(c(0,0),c(0,100),col="red")
    if(!to.file) browser()
  }
  ####################################################################################
  # deltas - fix everything and look at time
  ####################################################################################
  par(mfrow=c(5,2),mar=c(4,10,3,2))
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
    boxplot(y~x,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),main=type,xlab="log POD ratio",ylab="",ylim=c(-2,2))
    lines(c(0,0),c(0,100),col="red")
    if(!to.file) browser()
  }
  ####################################################################################
  # by condition
  ####################################################################################
  par(mfrow=c(1,1),mar=c(4,5,6,2))
  for(condition in unique(res$condition)) {
    temp = res[is.element(res$condition,condition),]
    temp = temp[order(temp$signature_pod_95),]
    ptemp = plier[is.element(plier$condition,condition),]
    rownames(ptemp) = ptemp$dtxsid

    plot(c(1,1),type="n",main=paste("By Condition\n",condition),xlim=c(1e-6,1000),ylim=c(0,44),cex.axis=1.2,cex.lab=1.2,xlab="POD (uM)",ylab="",log="x")
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
    plot(c(1,1),type="n",main=paste("By Type\n",type),xlim=c(1e-6,1000),ylim=c(0,44),cex.axis=1.2,cex.lab=1.2,xlab="POD (uM)",ylab="",log="x")
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
  boxplot(y~x,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),main="Range of PODs\nby Media and Time",xlab="log range",ylab="",ylim=c(0,4))

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
  boxplot(y~x,horizontal=T,las=1,par(cex.lab=1,cex.axis=1.0),main="Range of PODs\nby POD Type",xlab="log range",ylab="",ylim=c(0,4))


  if(!to.file) browser()
  if(to.file) dev.off()
}
