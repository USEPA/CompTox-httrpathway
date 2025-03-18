#--------------------------------------------------------------------------------------
#'
#' Build lane plots by chemical list and signature class, across the datasets
#'
#' @param to.file If TRUE, write plots to a file
#' @param dataset The data set to use
#' @param sigset THe signature set to use
#' @param method Scoring method
#' @param hccut Exclude rows with hitcall less than this value
#' @param plot.signature_min If TRUE, plot the minimum signature
#' @param bmd.mode percent or abs
#' @importFrom openxlsx read.xlsx write.xlsx
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics par points text lines
#' @importFrom stats lm
#' @export podLaneplot
#--------------------------------------------------------------------------------------
podLaneplot <- function(to.file=F,
                        dataset="DMEM_6hr_pilot_normal_pe_1",
                        sigset="pilot_large_all_100CMAP",
                        method="gsea",
                        hccut=0.9,
                        plot.signature_min=F,
                        bmd.mode="percent") {
  printCurrentFunction()
  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[order(chems$name),]
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)

  if(to.file) {
    fname <- paste0("../output/signature_pod/pod_laneplot_",dataset,"_",sigset,"_",method,"_",bmd.mode,".pdf")
    if(!plot.signature_min) fname <- paste0("../output/signature_pod/pod_laneplot_",dataset,"_",sigset,"_",method,"_",bmd.mode,"_no_signature_min.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,2))

  file <- "../toxcast/toxcast_pod.xlsx"
  toxcast <- read.xlsx(file)
  rownames(toxcast) <- toxcast$dtxsid

  file <- "../input/BMDExpress/BMDExpress_Pseudo1_ANOVA_Pathway_Results_Pilot_6h_DMEM.RDS"
  all_pathway_bmds <- readRDS(file) 
  ### all_pathway_bmds
  all_pathway_bmds <- all_pathway_bmds[is.element(all_pathway_bmds$type,"Real"),]

  file <- paste0("../output/signature_pod/signature_pod_",sigset,"_",dataset,"_",method,"_",hccut,".xlsx")
  print(file)
  pod.signature <- read.xlsx(file)
  if(bmd.mode=="abs") {
    pod.signature <- pod.signature[order(pod.signature$signature_pod_95.count),]
    #pod.signature = pod.signature[pod.signature$signature_pod_95.count<1000,]
  }
  else if(bmd.mode=="percent") pod.signature <- pod.signature[order(pod.signature$signature_pod_95),]

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
  plot(c(1,1),type="n",main=paste(dataset,":",method,"\n",bmd.mode,":",sigset),
       cex.main=1.0,cex.axis=1.2,cex.lab=1.2,xlab="log(BMD uM)",ylab="",yaxt="n",
       xlim=c(xmin,xmax),ylim=c(0,nchem+2),log="x",
       xaxp=c(1e-4,1e2,n=1))

  yval <- -0.5
  yval <- nchem+2.5
  xval <- 1.8e-7
  text(xval,yval,"5th Signature",pos=4,cex=1)
  points(xval,yval,pch=24,bg="black")
  xval <- xval*75
  if(plot.signature_min) {
    text(xval,yval,"Pathway min",pos=4,cex=1)
    points(xval,yval,pch=24,bg="blue")
    xval <- xval*75
  }
  points(xval,yval,pch=23,bg="red")
  text(xval,yval,"ToxCast",pos=4,cex=1)
  xval <- xval*20
  points(xval,yval,pch=23,bg="yellow")
  text(xval,yval,"BMDExpress",pos=4,cex=1)
  xval <- xval*80
  points(xval,yval,pch=24,bg="green")
  xval <- xval*2
  points(xval,yval,pch=25,bg="green")
  text(xval,yval,"ER agonist/antagonist",pos=4,cex=1)

  counter <- 0
  lines(c(xmin,xmax),c(counter+0.5,counter+0.5),col="gray")
  counter <- 1

  name.list <- c("dtxsid","casrn","name",
                 "toxcast_pod",
                 "signature_bmdl","signature_bmd","signature_bmdu","signature_pod_min",
                 "bmdexpress","er_agonist","er_antagonist")
  result <- as.data.frame(matrix(nrow=length(dtxsid.list),ncol=length(name.list)))
  names(result) <- name.list

  result$dtxsid <- chems$dtxsid
  result$casrn <- chems$casrn
  result$name <- chems$name
  rownames(result) <- result$dtxsid

  for(dtxsid in dtxsid.list) {
    name <- chems[is.element(chems$dtxsid,dtxsid),"name"]
    lines(c(xmin,xmax),c(counter+0.5,counter+0.5),col="gray")

    signature_pod_min <- pod.signature[dtxsid,"signature_pod_min"]
    signature_pod_min.lci <- pod.signature[dtxsid,"signature_pod_min.lci"]
    signature_pod_min.uci <- pod.signature[dtxsid,"signature_pod_min.uci"]
    if(bmd.mode=="abs") {
      signature_pod_95 <- pod.signature[dtxsid,"signature_pod_95.count"]
      signature_pod_95.lci <- pod.signature[dtxsid,"signature_pod_95.count.lci"]
      signature_pod_95.uci <- pod.signature[dtxsid,"signature_pod_95.count.uci"]
    }
    else if(bmd.mode=="percent") {
      signature_pod_95 <- pod.signature[dtxsid,"signature_pod_95"]
      signature_pod_95.lci <- pod.signature[dtxsid,"signature_pod_95.lci"]
      signature_pod_95.uci <- pod.signature[dtxsid,"signature_pod_95.uci"]
    }

    toxcast_pod_05 <- toxcast[dtxsid,"pod_uM"]

    bmds_pod = 10000
    if(is.element(name,all_pathway_bmds$chem_name))
      bmds_pod <- min(all_pathway_bmds[is.element(all_pathway_bmds$chem_name,name),"BMD"],na.rm=T)

    if(is.na(signature_pod_min)) signature_pod_min <- 1000
    if(is.na(signature_pod_min.lci)) signature_pod_min.lci <- 0.001
    if(is.na(signature_pod_min.uci)) signature_pod_min.uci <- 1000
    if(is.na(signature_pod_95)) signature_pod_95 <- 1000
    if(is.na(signature_pod_95.lci)) signature_pod_95.lci <- 0.001
    if(is.na(signature_pod_95.uci)) signature_pod_95.uci <- 1000
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

    points(toxcast_pod_05,counter,pch=23,bg="red")

    if(!is.na(bmds_pod)) points(bmds_pod[1],counter,pch=23,bg="yellow")

    if(plot.signature_min) {
      if(signature_pod_min.lci<1e-4) signature_pod_min.lci <- 1e-4
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

    shift <- "="
    if(toxcast_pod_05<signature_pod_95.lci) shift <- "-"
    if(toxcast_pod_05>signature_pod_95.uci) shift <- "+"
    sname <- name
    col <- "black"
    if(shift=="=") col <- "red"
    if(shift=="+") col <- "royalblue"
    text(1e-4,counter,sname,pos=2,cex=0.8,col=col)

    result[dtxsid,"toxcast_pod"] <- toxcast_pod_05
    result[dtxsid,"signature_bmdl"] <- signature_pod_95.lci
    result[dtxsid,"signature_bmd"] <- signature_pod_95
    result[dtxsid,"signature_bmdu"] <- signature_pod_95.uci
    result[dtxsid,"signature_pod_min"] <- signature_pod_min
    result[dtxsid,"bmdexpress"] <- bmds_pod

    counter <- counter+1
  }
  yval <- -0.5
  xval <- 0.5e-7

  x <- log10(result$signature_bmd)
  y <- log10(result$toxcast_pod)
  res <- lm(y~x)
  sr <- summary(res)
  r2.1 <- sr$adj.r.squared
  rmse.1 <- sr$sigma

  x <- log10(result$signature_bmd)
  y <- log10(result$bmdexpress)
  y[y==Inf] <- 1000
  res <- lm(y~x)
  sr <- summary(res)
  r2.2 <- sr$adj.r.squared
  rmse.2 <- sr$sigma

  footer <- paste0("black,red,blue: ToxCast POD <,in,> Pathway POD range; R2(ToxCast, BMDExpress)=",
                   format(r2.1,digits=2),":",format(r2.2,digits=2))

  text(xval,yval,footer,pos=4,cex=0.8)
  if(!to.file) browser()
  else dev.off()

  if(!plot.signature_min) file <- paste0("../output/signature_pod/pod_laneplot_",dataset,"_",sigset,"_",method,"_",bmd.mode,".xlsx")
  else file <- paste0("../output/signature_pod/pod_laneplot_",dataset,"_",sigset,"_",method,"_",bmd.mode,"_no_signature_min.xlsx")
  write.xlsx(result,file)
}

