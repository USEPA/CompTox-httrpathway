#--------------------------------------------------------------------------------------
#'
#' Build lane plots by chemical list and pathway class, across the datasets
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
podLaneplot <- function(to.file=F,
                        dataset="DMEM_6hr_pilot_normal_pe_0",
                        pathset="PathwaySet_20191107",
                        method="fc") {
  printCurrentFunction()
  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[order(chems$name),]
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)

  if(to.file) {
    fname <- paste0("../output/pod_laneplot/pod_laneplot_",dataset,"_",pathset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(4,4,2,2))

  file <- paste0("../output/accumplots/NULLACCUMPLOT_",pathset,"_",dataset,"_0.05_conthits.xlsx")
  print(file)
  pod.accum <- read.xlsx(file)
  rownames(pod.accum) <- pod.accum$dtxsid

  file <- paste0("../output/pathway_pod/pathway_pod_",pathset,"_",dataset,"_",method,".xlsx")
  print(file)
  pod.pathway <- read.xlsx(file)
  pod.pathway <- pod.pathway[order(pod.pathway$pathway_pod_95),]
  rownames(pod.pathway) <- pod.pathway$dtxsid
  dtxsid.list <- pod.pathway$dtxsid

  file <- "../input/ER/ER_model_data 2019-11-01.xlsx"
  erdata <- read.xlsx(file)
  rownames(erdata) <- erdata$dtxsid

  xmin <- 1e-7
  xmax <- 1e3
  plot(c(1,1),type="n",main="POD Summary",cex.axis=1.2,cex.lab=1.2,xlab="log(BMD10 uM)",ylab="",yaxt="n",
       xlim=c(xmin,xmax),ylim=c(0,nchem+1),log="x",
       xaxp=c(1e-4,1e2,n=1))
  points(1e-6,-0.5,pch=21,bg="red")
  text(1e-6,-0.5,"pathway acc",pos=4,cex=1)
  points(1e-4,-0.5,pch=24,bg="orange")
  text(1e-4,-0.5,"pathway LCI",pos=4,cex=1)
  points(1e-2,-0.5,pch=24,bg="blue")
  text(1e-2,-0.5,"pathway min",pos=4,cex=1)
  points(1e-0,-0.5,pch=22,bg="green")
  text(1e-0,-0.5,"ER agonist",pos=4,cex=1)

  counter <- 0
  lines(c(xmin,xmax),c(counter+0.5,counter+0.5),col="gray")
  counter <- 1
  for(dtxsid in dtxsid.list) {
    name <- chems[is.element(chems$dtxsid,dtxsid),"name"]
    lines(c(xmin,xmax),c(counter+0.5,counter+0.5),col="gray")
    text(1e-4,counter,name,pos=2,cex=0.8)
    accum.bmd <- pod.accum[dtxsid,"bmd"]
    pathway_pod_min <- pod.pathway[dtxsid,"pathway_pod_min"]
    pathway_pod_95 <- pod.pathway[dtxsid,"pathway_pod_95"]

    if(is.na(accum.bmd)) accum.bmd <- 1000
    if(is.na(pathway_pod_min)) pathway_pod_min <- 1100
    if(is.na(pathway_pod_95)) pathway_pod_95 <- 1200
    points(accum.bmd,counter,pch=21,bg="red")
    points(pathway_pod_min,counter,pch=24,bg="blue")
    points(pathway_pod_95,counter,pch=24,bg="orange")

    if(is.element(dtxsid,erdata$dtxsid)) {
      erpod <- 10**(erdata[dtxsid,"AC50meanFromAUCcal"])
      if(erpod<100) points(erpod,counter,pch=22,bg="green")
    }

    counter <- counter+1
  }
  if(!to.file) browser()
  else dev.off()
}

