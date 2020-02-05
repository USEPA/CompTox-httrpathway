#--------------------------------------------------------------------------------------
#'
#' Build lane plots by chemical list and signature class, across the datasets
#' @param to.file If TRUE, write plots to a file
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#--------------------------------------------------------------------------------------
signatureChemicalLanePlotAcrossDatasets <- function(to.file=F,
                                                  chemical.target="ER",
                                                  signature.super_class="estrogen",
                                                  pathset="PathwaySet_20191107",
                                                  method = "fc") {
  printCurrentFunction()
  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[is.element(chems$target_key,chemical.target),]
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)

  file <- "../input/processed_signature_data/signature_catalog 2019-11-07.xlsx"
  catalog <- read.xlsx(file)
  catalog <- catalog[catalog$useme==1,]
  catalog <- catalog[is.element(catalog$super_class,signature.super_class),]
  signature.list <- catalog$signature
  signature.list <- sort(signature.list)

  dataset.list = c(
    "DMEM_6hr_pilot_normal_pe_0",
    "DMEM_12hr_pilot_normal_pe_0",
    "DMEM_24hr_pilot_normal_pe_0",

    "DMEM_6hr_pilot_none_pe_0",
    "DMEM_12hr_pilot_none_pe_0",
    "DMEM_24hr_pilot_none_pe_0",

    "DMEM_6hr_pilot_normal_pe_1",
    "DMEM_12hr_pilot_normal_pe_1",
    "DMEM_24hr_pilot_normal_pe_1",

    "DMEM_6hr_pilot_none_pe_1",
    "DMEM_12hr_pilot_none_pe_1",
    "DMEM_24hr_pilot_none_pe_1",

    "PRF_6hr_pilot_normal_pe_0",
    "PRF_12hr_pilot_normal_pe_0",
    "PRF_24hr_pilot_normal_pe_0",

    "PRF_6hr_pilot_none_pe_0",
    "PRF_12hr_pilot_none_pe_0",
    "PRF_24hr_pilot_none_pe_0",

    "PRF_6hr_pilot_normal_pe_1",
    "PRF_12hr_pilot_normal_pe_1",
    "PRF_24hr_pilot_normal_pe_1",

    "PRF_6hr_pilot_none_pe_1",
    "PRF_12hr_pilot_none_pe_1",
    "PRF_24hr_pilot_none_pe_1"
  )
  if(to.file) {
    fname <- paste0("../output/signature_conc_resp_laneplots/signatureChemicalLanePlotAcrossDatasets_",chemical.target,"_",signature.super_class,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,1),mar=c(4,4,2,2))

  dataset.list <- dataset.list[1:12]
  #dataset.list <- dataset.list[1]

  mat <- NULL

  #signature.list <- signature.list[1:3]
  for(dataset in dataset.list) {
    file <- paste0("../output/signature_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    temp <- PATHWAY_CR
    temp <- temp[is.element(temp$dtxsid,dtxsid.list),]
    temp <- temp[is.element(temp$signature,signature.list),]
    temp$dataset <- dataset
    mat <- rbind(mat,temp)
  }

  xmin <- -6
  xmax <- 3
  for(signature in signature.list) {
    cat(signature,"\n")
    plot(c(1,1),type="n",main=signature,cex.axis=1.2,cex.lab=1.2,xlab="log(BMD uM)",ylab="",yaxt="n",
         xlim=c(xmin,xmax),ylim=c(0,nchem+1))
    points(-5,0,pch=21,bg="red")
    text(-5,0,"normal",pos=4,cex=1.2)
    points(-4,0,pch=21,bg="cyan")
    text(-4,0,"none",pos=4,cex=1.2)
    points(-3,0,pch=21,bg="white")
    text(-3,0,"pe=0",pos=4,cex=1.2)
    points(-2,0,pch=24,bg="white")
    text(-2,0,"pe=1",pos=4,cex=1.2)

    temp1 <- mat[is.element(mat$signature,signature),]
    counter <- 0
    lines(c(xmin,xmax),c(counter+0.5,counter+0.5),col="gray")
    counter <- 1
    for(dtxsid in dtxsid.list) {
      name <- chems[is.element(chems$dtxsid,dtxsid),"name"]
      lines(c(xmin,xmax),c(counter+0.5,counter+0.5),col="gray")
      text(-4,counter,name,pos=2,cex=1.2)
      temp2 <- temp1[is.element(temp1$dtxsid,dtxsid),]
      x <- temp2$bmd
      y <- temp2$hitcall
      x[is.na(x)] <- 1000
      x[y<0.5] <- 1000
      temp2$bmd <- x
      for(dataset in dataset.list) {
        offset <- 0
        if(contains(dataset,"6h")) offset <- -0.3
        if(contains(dataset,"24h")) offset <- 0.3
        color <- "gray"
        if(contains(dataset,"normal")) color <- "red"
        if(contains(dataset,"none")) color <- "cyan"
        pch <- 21
        if(contains(dataset,"pe_1")) pch <- 24
        bmd <- temp2[is.element(temp2$dataset,dataset),"bmd"]
        #cat(dataset,pch,bmd,"\n")
        if(bmd==1000) bmd <- 1000+rnorm(1,0,50)
        bmd <- log10(bmd)
        points(bmd,counter+offset,pch=pch,cex=1.5,bg=color)
      }
      counter <- counter+1

    }

    if(!to.file) browser()
  }
  if(to.file) dev.off()
}

