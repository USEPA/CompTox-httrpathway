#--------------------------------------------------------------------------------------
#'
#' Build lane plots by chemical list and signature class, across the datasets
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
podLaneplot <- function(to.file=F,
                        dataset,
                        sigset,
                        method,
                        plot.signature_min=F) {
  printCurrentFunction()
  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[order(chems$name),]
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)

  if(to.file) {
    fname <- paste0("../output/signature_pod/pod_laneplot_",dataset,"_",sigset,"_",method,".pdf")
    if(!plot.signature_min) fname <- paste0("../output/signature_pod/pod_laneplot_",dataset,"_",sigset,"_",method,"_no_signature_min.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,2))

  file <- "../toxcast/toxcast_pod.xlsx"
  print(file)
  toxcast <- read.xlsx(file)
  rownames(toxcast) <- toxcast$dtxsid

  #file <- paste0("../output/accumplots/NULLACCUMPLOT_",sigset,"_",dataset,"_0.05_conthits.xlsx")
  #print(file)
  #pod.accum <- read.xlsx(file)
  #rownames(pod.accum) <- pod.accum$dtxsid

  file <- "../input/BMDExpress/BMDExpress_Pathway_Results_Pilot_6h_DMEM.RData"
  load(file=file)
  ### min_pathway_bmds

  file <- paste0("../output/signature_pod/signature_pod_",sigset,"_",dataset,"_",method,".xlsx")
  print(file)
  pod.signature <- read.xlsx(file)
  pod.signature <- pod.signature[order(pod.signature$signature_pod_95),]
  rownames(pod.signature) <- pod.signature$dtxsid
  dtxsid.list <- pod.signature$dtxsid

  file <- "../input/ER/lcia_ar_er_data.xlsx"
  erdata <- read.xlsx(file)
  rownames(erdata) <- erdata$dtxsid

  file <- "../input/cytotoxicity summary wide phase12.xlsx"
  cytotox <- read.xlsx(file)
  rownames(cytotox) <- cytotox$dtxsid

  xmin <- 1e-7
  xmax <- 1e3
  plot(c(1,1),type="n",main=paste(dataset,":",method,"\n",sigset),cex.main=1.0,cex.axis=1.2,cex.lab=1.2,xlab="log(BMD uM)",ylab="",yaxt="n",
       xlim=c(xmin,xmax),ylim=c(0,nchem+2),log="x",
       xaxp=c(1e-4,1e2,n=1))

  yval <- -0.5
  yval <- nchem+2.5
  xval <- 1.8e-7
  text(xval,yval,"5th Signature",pos=4,cex=1)
  points(xval,yval,pch=24,bg="black")
  xval <- xval*100
  if(plot.signature_min) {
    text(xval,yval,"Pathway min",pos=4,cex=1)
    points(xval,yval,pch=24,bg="blue")
    xval <- xval*100
  }
  #points(xval,yval,pch=21,bg="orange")
  #text(xval,yval,"Pathway acc",pos=4,cex=1)
  #xval <- xval*100
  points(xval,yval,pch=23,bg="red")
  text(xval,yval,"ToxCast",pos=4,cex=1)
  xval <- xval*40
  points(xval,yval,pch=23,bg="yellow")
  text(xval,yval,"BMDExpress",pos=4,cex=1)
  xval <- xval*100
  points(xval,yval,pch=24,bg="green")
  text(xval,yval,"ER agonist",pos=4,cex=1)
  xval <- xval*50
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
                 "signature_bmdl","signature_bmd","signature_bmdu",
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
    #accum.bmd <- pod.accum[dtxsid,"bmd"]

    signature_pod_min <- pod.signature[dtxsid,"signature_pod_min"]
    signature_pod_min.lci <- pod.signature[dtxsid,"signature_pod_min.lci"]
    signature_pod_min.uci <- pod.signature[dtxsid,"signature_pod_min.uci"]

    signature_pod_95 <- pod.signature[dtxsid,"signature_pod_95"]
    signature_pod_95.lci <- pod.signature[dtxsid,"signature_pod_95.lci"]
    signature_pod_95.uci <- pod.signature[dtxsid,"signature_pod_95.uci"]

    toxcast_pod_05 <- toxcast[dtxsid,"pod_uM"]

    bmds_pod <- min(min_pathway_bmds[is.element(min_pathway_bmds$chem_name,name),"BMD"])

    if(is.na(signature_pod_min)) signature_pod_min <- 1000
    if(is.na(signature_pod_min.lci)) signature_pod_min.lci <- 0.001
    if(is.na(signature_pod_min.uci)) signature_pod_min.uci <- 1000
    if(is.na(signature_pod_95)) signature_pod_95 <- 1000
    if(is.na(signature_pod_95.lci)) signature_pod_95.lci <- 0.001
    if(is.na(signature_pod_95.uci)) signature_pod_95.uci <- 1000
    #if(is.na(accum.bmd)) accum.bmd <- 1000
    if(is.na(toxcast_pod_05)) toxcast_pod_05 <- 1000

    if(signature_pod_95.uci==1000 && signature_pod_95.lci==0.001) signature_pod_95.lci <- 1000
    if(signature_pod_min.uci==1000 && signature_pod_min.lci==0.001) signature_pod_min.lci <- 1000

    delta <- 0.5
    xm <- log10(signature_pod_95.lci)
    x0 <- log10(signature_pod_95)
    xp <- log10(signature_pod_95.uci)
    if(x0-xm < delta) xm <- x0-delta
    if(xp-x0 < delta) xp <- x0+delta
    signature_pod_95.lci <- 10**xm
    signature_pod_95.uci <- 10**xp

    #points(accum.bmd,counter,pch=21,bg="orange")
    points(toxcast_pod_05,counter,pch=23,bg="red")

    if(!is.na(bmds_pod)) points(bmds_pod[1],counter,pch=23,bg="yellow")

    if(plot.signature_min) {
      if(signature_pod_min.lci<1e-4) signature_pod_min.lci <- 1e-4
      lines(c(signature_pod_min.lci,signature_pod_min.uci),c(counter-0.25,counter-0.25),col="blue",lwd=2)
      points(signature_pod_min,counter-0.25,pch=24,bg="blue")
    }
    lines(c(signature_pod_95.lci,signature_pod_95.uci),c(counter+0.25,counter+0.25),col="black",lwd=2)
    points(signature_pod_95,counter+0.25,pch=24,bg="black")

    if(is.element(dtxsid,erdata$dtxsid)) {
      auc.agon <- erdata[dtxsid,"ER.agonist.auc"]
      auc.anta <- erdata[dtxsid,"ER.antagonist.auc"]
      if(is.na(auc.agon)) auc.agon <- 0
      if(is.na(auc.anta)) auc.anta <- 0
      if(auc.agon>auc.anta) {
        erpod <- 10**(erdata[dtxsid,"ER.agonist.meanlogac50.mean"])
        if(erpod<100 && name!="Reserpine") {
          points(erpod,counter,pch=24,bg="green")
          result[dtxsid,"er_agonist"] <- erpod
        }
      }
      else {
        erpod <- 10**(erdata[dtxsid,"ER.antagonist.meanlogac50.mean"])
        if(erpod<100 && name!="Reserpine") {
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
    #if(cytomin < signature_pod_95 ) scyto <- "X"
    text(2e-7,counter,scyto,pos=2,cex=0.8)

    shift <- "="
    if(toxcast_pod_05<signature_pod_95.lci) shift <- "-"
    if(toxcast_pod_05>signature_pod_95.uci) shift <- "+"
    sname <- paste0(scyto," ",name," (",shift,")")
    sname <- name
    col <- "black"
    if(shift=="=") col <- "red"
    if(shift=="+") col <- "royalblue"
    text(1e-4,counter,sname,pos=2,cex=0.8,col=col)

    result[dtxsid,"toxcast_pod"] <- toxcast_pod_05
    result[dtxsid,"signature_bmdl"] <- signature_pod_95.lci
    result[dtxsid,"signature_bmd"] <- signature_pod_95
    result[dtxsid,"signature_bmdu"] <- signature_pod_95.uci
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

  file <- paste0("../output/signature_pod/pod_laneplot_",dataset,"_",sigset,"_",method,"_no_signature_min.xlsx")
  write.xlsx(result,file)

}

