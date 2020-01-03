#--------------------------------------------------------------------------------------
#'
#' Build lane plots by chemical list and pathway class, across the datasets
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
podLaneplot <- function(to.file=F,
                        dataset="DMEM_6hr_pilot_normal_pe_1",
                        pathset="PathwaySet_20191107",
                        method="gsva",
                        plot.pathway_min=F) {
  printCurrentFunction()
  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[order(chems$name),]
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)

  if(to.file) {
    fname <- paste0("../output/pod_laneplot/pod_laneplot_",dataset,"_",pathset,"_",method,".pdf")
    if(!plot.pathway_min) fname <- paste0("../output/pod_laneplot/pod_laneplot_",dataset,"_",pathset,"_",method,"_no_pathway_min.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,2))

  file <- "../toxcast/toxcast_pod.xlsx"
  print(file)
  toxcast <- read.xlsx(file)
  rownames(toxcast) <- toxcast$dtxsid

  file <- paste0("../output/accumplots/NULLACCUMPLOT_",pathset,"_",dataset,"_0.05_conthits.xlsx")
  print(file)
  pod.accum <- read.xlsx(file)
  rownames(pod.accum) <- pod.accum$dtxsid

  file <- "../input/BMDExpress/BMDExpress_Pathway_Results_Pilot_6h_DMEM.RData"
  load(file=file)
  ### min_pathway_bmds

  file <- paste0("../output/pathway_pod/pathway_pod_",pathset,"_",dataset,"_",method,".xlsx")
  print(file)
  pod.pathway <- read.xlsx(file)
  pod.pathway <- pod.pathway[order(pod.pathway$pathway_pod_95),]
  rownames(pod.pathway) <- pod.pathway$dtxsid
  dtxsid.list <- pod.pathway$dtxsid

  file <- "../input/ER/lcia_ar_er_data.xlsx"
  erdata <- read.xlsx(file)
  rownames(erdata) <- erdata$dtxsid

  file <- "../input/cytotoxicity summary wide phase12.xlsx"
  cytotox <- read.xlsx(file)
  rownames(cytotox) <- cytotox$dtxsid

  xmin <- 1e-7
  xmax <- 1e3
  plot(c(1,1),type="n",main=paste(dataset,":",method,"\n",pathset),cex.main=1.0,cex.axis=1.2,cex.lab=1.2,xlab="log(BMD10 uM)",ylab="",yaxt="n",
       xlim=c(xmin,xmax),ylim=c(0,nchem+2),log="x",
       xaxp=c(1e-4,1e2,n=1))

  yval <- -0.5
  yval <- nchem+2.5
  xval <- 1.8e-7
  text(xval,yval,"Pathway LCI",pos=4,cex=1)
  points(xval,yval,pch=24,bg="black")
  xval <- xval*100
  if(plot.pathway_min) {
    text(xval,yval,"Pathway min",pos=4,cex=1)
    points(xval,yval,pch=24,bg="blue")
    xval <- xval*100
  }
  #points(xval,yval,pch=21,bg="orange")
  #text(xval,yval,"Pathway acc",pos=4,cex=1)
  #xval <- xval*100
  points(xval,yval,pch=23,bg="red")
  text(xval,yval,"ToxCast",pos=4,cex=1)
  xval <- xval*100
  points(xval,yval,pch=23,bg="yellow")
  text(xval,yval,"BMDExpress",pos=4,cex=1)
  xval <- xval*100
  points(xval,yval,pch=24,bg="green")
  text(xval,yval,"ER agonist",pos=4,cex=1)
  xval <- xval*100
  points(xval,yval,pch=25,bg="green")
  text(xval,yval,"ER antagonist",pos=4,cex=1)

  xval <- 1e-7
  yval <- nchem+1.25
  #text(xval,yval,"Cytotoxicity: ",pos=4); xval<- xval*100
  #text(xval,yval,"BLA",pos=4)
  #lines(c(xval,xval),c(yval-0.5,yval+0.5),lwd=2,col="red"); xval <- xval*10
  #text(xval,yval,"LUC",pos=4)
  #lines(c(xval,xval),c(yval-0.5,yval+0.5),lwd=2,col="orange"); xval <- xval*10
  #text(xval,yval,"SRB",pos=4)
  #lines(c(xval,xval),c(yval-0.5,yval+0.5),lwd=2,col="cyan"); xval <- xval*10
  #text(xval,yval,"FLO",pos=4)
  #lines(c(xval,xval),c(yval-0.5,yval+0.5),lwd=2,col="gray"); xval <- xval*10
  #text(xval,yval,"Other",pos=4)
  #lines(c(xval,xval),c(yval-0.5,yval+0.5),lwd=2,col="blue"); xval <- xval*10
  counter <- 0
  lines(c(xmin,xmax),c(counter+0.5,counter+0.5),col="gray")
  counter <- 1

  name.list <- c("dtxsid","casrn","name",
                 "toxcast_pod",
                 "pathway_bmdl","pathway_bmd","pathway_bmdu",
                 "bmdexpress","er_agonist","er_antagonist",
                 "cytotox_bla","cytotox_luc","cytotox_srb","cytotox_other")
  result <- as.data.frame(matrix(nrow=length(dtxsid.list),ncol=length(name.list)))
  names(result) <- name.list

  result$dtxsid <- chems$dtxsid
  result$casrn <- chems$casrn
  result$name <- chems$name
  rownames(result) <- result$dtxsid

  for(dtxsid in dtxsid.list) {
    name <- chems[is.element(chems$dtxsid,dtxsid),"name"]
    lines(c(xmin,xmax),c(counter+0.5,counter+0.5),col="gray")
    accum.bmd <- pod.accum[dtxsid,"bmd"]

    pathway_pod_min <- pod.pathway[dtxsid,"pathway_pod_min"]
    pathway_pod_min.lci <- pod.pathway[dtxsid,"pathway_pod_min.lci"]
    pathway_pod_min.uci <- pod.pathway[dtxsid,"pathway_pod_min.uci"]

    pathway_pod_95 <- pod.pathway[dtxsid,"pathway_pod_95"]
    pathway_pod_95.lci <- pod.pathway[dtxsid,"pathway_pod_95.lci"]
    pathway_pod_95.uci <- pod.pathway[dtxsid,"pathway_pod_95.uci"]

    toxcast_pod_05 <- toxcast[dtxsid,"pod_uM"]

    bmds_pod <- min(min_pathway_bmds[is.element(min_pathway_bmds$chem_name,name),"BMD"])

    if(is.na(pathway_pod_min)) pathway_pod_min <- 1000
    if(is.na(pathway_pod_min.lci)) pathway_pod_min.lci <- 0.001
    if(is.na(pathway_pod_min.uci)) pathway_pod_min.uci <- 1000
    if(is.na(pathway_pod_95)) pathway_pod_95 <- 1000
    if(is.na(pathway_pod_95.lci)) pathway_pod_95.lci <- 0.001
    if(is.na(pathway_pod_95.uci)) pathway_pod_95.uci <- 1000
    if(is.na(accum.bmd)) accum.bmd <- 1000
    if(is.na(toxcast_pod_05)) toxcast_pod_05 <- 1000

    if(pathway_pod_95.uci==1000 && pathway_pod_95.lci==0.001) pathway_pod_95.lci <- 1000
    if(pathway_pod_min.uci==1000 && pathway_pod_min.lci==0.001) pathway_pod_min.lci <- 1000

    delta <- 0.5
    xm <- log10(pathway_pod_95.lci)
    x0 <- log10(pathway_pod_95)
    xp <- log10(pathway_pod_95.uci)
    if(x0-xm < delta) xm <- x0-delta
    if(xp-x0 < delta) xp <- x0+delta
    pathway_pod_95.lci <- 10**xm
    pathway_pod_95.uci <- 10**xp

    #points(accum.bmd,counter,pch=21,bg="orange")
    points(toxcast_pod_05,counter,pch=23,bg="red")

    if(!is.na(bmds_pod)) points(bmds_pod[1],counter,pch=23,bg="yellow")

    if(plot.pathway_min) {
      if(pathway_pod_min.lci<1e-4) pathway_pod_min.lci <- 1e-4
      lines(c(pathway_pod_min.lci,pathway_pod_min.uci),c(counter-0.25,counter-0.25),col="blue",lwd=2)
      points(pathway_pod_min,counter-0.25,pch=24,bg="blue")
    }
    lines(c(pathway_pod_95.lci,pathway_pod_95.uci),c(counter+0.25,counter+0.25),col="black",lwd=2)
    points(pathway_pod_95,counter+0.25,pch=24,bg="black")

    if(is.element(dtxsid,erdata$dtxsid)) {
      auc.agon <- erdata[dtxsid,"ER.agonist.auc"]
      auc.anta <- erdata[dtxsid,"ER.antagonist.auc"]
      if(is.na(auc.agon)) auc.agon <- 0
      if(is.na(auc.anta)) auc.anta <- 0
      if(auc.agon>auc.anta) {
        erpod <- 10**(erdata[dtxsid,"ER.agonist.meanlogac50.mean"])
        if(erpod<100) {
          points(erpod,counter,pch=24,bg="green")
          result[dtxsid,"er_agonist"] <- erpod
        }
      }
      else {
        erpod <- 10**(erdata[dtxsid,"ER.antagonist.meanlogac50.mean"])
        if(erpod<100) {
          points(erpod,counter,pch=25,bg="green")
          result[dtxsid,"er_antagonist"] <- erpod
        }
      }
    }

    luc <- 10**cytotox[dtxsid,"LUC"]
    bla <- 10**cytotox[dtxsid,"BLA"]
    srb <- 10**cytotox[dtxsid,"SRB"]
    flo <- 10**cytotox[dtxsid,"FLO"]
    other <- 10**cytotox[dtxsid,"other"]
    #if(!is.na(luc)) lines(c(luc,luc),c(counter-0.5,counter+0.5),col="orange",lwd=2)
    #if(!is.na(bla)) lines(c(bla,bla),c(counter-0.5,counter+0.5),col="red",lwd=2)
    #if(!is.na(srb)) lines(c(srb,srb),c(counter-0.5,counter+0.5),col="cyan",lwd=2)
    #if(!is.na(flo)) lines(c(flo,flo),c(counter-0.5,counter+0.5),col="gray",lwd=2)
    #if(!is.na(other)) lines(c(other,other),c(counter-0.5,counter+0.5),col="blue",lwd=2)
    cytomin <- min(luc,bla,srb,other)
    scyto <- ""
    #if(cytomin < pathway_pod_95 ) scyto <- "X"
    text(2e-7,counter,scyto,pos=2,cex=0.8)

    shift <- "="
    if(toxcast_pod_05<pathway_pod_95.lci) shift <- "-"
    if(toxcast_pod_05>pathway_pod_95.uci) shift <- "+"
    sname <- paste0(scyto," ",name," (",shift,")")
    sname <- name
    col <- "black"
    if(shift=="=") col <- "red"
    if(shift=="+") col <- "royalblue"
    text(1e-4,counter,sname,pos=2,cex=0.8,col=col)

    result[dtxsid,"toxcast_pod"] <- toxcast_pod_05
    result[dtxsid,"pathway_bmdl"] <- pathway_pod_95.lci
    result[dtxsid,"pathway_bmd"] <- pathway_pod_95
    result[dtxsid,"pathway_bmdu"] <- pathway_pod_95.uci
    result[dtxsid,"bmdexpress"] <- bmds_pod
    result[dtxsid,"cytotox_bla"] <- bla
    result[dtxsid,"cytotox_luc"] <- luc
    result[dtxsid,"cytotox_srb"] <- srb
    result[dtxsid,"cytotox_other"] <- other

    counter <- counter+1
  }
  yval <- -0.5
  xval <- 0.5e-7
  #text(xval,yval,"black,red,blue: ToxCast POD <,in,> Pathway POD range;  X: Cytotoxicity and Pathway POD range overlap",pos=4,cex=0.8)
  text(xval,yval,"black,red,blue: ToxCast POD <,in,> Pathway POD range",pos=4,cex=0.8)
  if(!to.file) browser()
  else dev.off()

  file <- paste0("../output/pod_laneplot/pod_laneplot_",dataset,"_",pathset,"_",method,"_no_pathway_min.xlsx")
  write.xlsx(result,file)

}

