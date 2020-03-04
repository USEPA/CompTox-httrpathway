#--------------------------------------------------------------------------------------
#'
#' Build summary plots by signature class
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
signatureClassSummaryPlot <- function(to.file=F,dataset="DMEM_6hr_pilot_normal_pe_1",
                                    sigset="pilot_large_all_CMAP",
                                    method = "mygsea",
                                    sigcatalog="signatureDB_master_catalog 2020-01-31",
                                    hitcall.threshold=0.5) {
  printCurrentFunction()

  file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  load(file=file)
  mat <- SIGNATURE_CR

  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  achems <- read.xlsx(file)
  rownames(achems) <- achems$dtxsid

  chems <- unique(mat[,c("dtxsid","casrn","name")])
  chems <- chems[order(chems$name),]
  mat <- mat[mat$hitcall>=hitcall.threshold,]
  mat <- mat[!is.na(mat$bmd),]

  if(to.file) {
    fname <- paste0("../output/signature_class_summary_plots/signatureClassSummaryPlot ",dataset,"_",sigset,"_",method,"_0.05_conthits.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  annotations <- read.xlsx(file)

  delta <- 0.25
  grid <- seq(from=-5,to=4,by=delta)

  x <- log10(mat$bmd)
  ix <- x
  ix[] <- 0
  for(i in 1:length(grid)) {
    ix[x>grid[i]] <- i
  }
  mat$logbmd <- x
  mat$bmdindex <- ix
  parent.list <- unique(mat$signature)
  colors <- mat$signature
  colors[] <- "white"
  for(parent in parent.list) {
    color <- annotations[is.element(annotations$parent,parent),"color"][1]
    #cat(parent,color,"\n")
    colors[is.element(mat$signature,parent)] <- color
  }

  mat$color <- colors
  mask <- mat$bmd
  mask[] <- 0
  mask1 <- mask
  mask1[mat$fit_method=="gnls"] <- 1
  mask2 <- mask
  mask2[mat$bmd< 0.1] <- 1
  mask <- 1- mask1*mask2
  mat <- mat[mask==1,]

  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chem.annotations <- read.xlsx(file)
  rownames(chem.annotations) <- chem.annotations$dtxsid
  npath <- length(unique(mat$signature))
  for(i in 1:nrow(chems)) {
    dtxsid <- chems[i,"dtxsid"]
    name <- chems[i,"name"]
    #dtxsid <- "DTXSID4022369"
    name <- chems[is.element(chems$dtxsid,dtxsid),"name"]
    ymin <- 0
    ymax <- 60

   # if(is.element(name,c("Cyproterone acetate","Maneb","Ziram"))) ymax <- 100
    plot(c(0,0),main=name,cex.axis=1.2,cex.lab=1.2,type="n",
         xlim=c(-4,2),ylim=c(0,ymax),xlab="log(bmd uM)",ylab="Pathway Efficacy")

    color <- achems[is.element(achems$dtxsid,dtxsid),"color"]
    rect(1.5,ymax-5,2,ymax,col=color,border="black")

    temp <- mat[is.element(mat$dtxsid,dtxsid),]
    cat("initial signatures",nrow(temp),"\n")
    temp <- mat[is.element(mat$dtxsid,dtxsid),c("signature","bmdindex","top","color","super_target")]

    target_class <- achems[dtxsid,"target_key"]
    target.list <- str_split(target_class,"\\|")[[1]]
    #color.list <- temp$color
    #color.list[!is.element(color.list,c("red","black"))] <- "gray"
    #color.list[is.element(temp$super_target,target.list)] <- "green"
    #temp$color <- color.list
    #browser()

    temp <- temp[order(temp$color),]
    #temp <- temp[!is.element(temp$color,"gray"),]
    cat(name,nrow(temp),"\n")
    itop <- grid
    itop[] <- 0
    for(j in 1:nrow(temp)) {
      index <- temp[j,"bmdindex"]
      x <- grid[index]
      y <- itop[index]
      top <- abs(temp[j,"top"])
      top <- 0.1
      color <- temp[j,"color"]
      rect(x,y,x+delta,y+top,col=color,border=color)
      itop[index] <- itop[index]+top

    }
    cclass <- chem.annotations[dtxsid,"target_annotation"]
    text(-4,ymax*0.95,cclass,pos=4,cex=0.95)
    text(-4,ymax*0.85,paste0("Pathways with hitcall>",hitcall.threshold,": ",nrow(temp)," / ",npath),pos=4,cex=0.95)

    dd <- (ymax-ymin)/8
    yval <- ymax * 0.8#ymin
    if(!to.file) browser()
  }
  if(to.file) dev.off()
}

