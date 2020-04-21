#--------------------------------------------------------------------------------------
#'
#' Look at the effect of combining the up and down version of the CMAP genesets
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
CMAP.up.dn <- function(to.file=F,
                       pathset.split="PathwaySet_20191107",
                       pathset.combined="PathwaySet_20200108",
                       dataset="DMEM_6hr_pilot_normal_pe_1",
                       method = "mygsea") {
  printCurrentFunction()
  if(to.file) {
    fname <- paste0("../output/pathway_conc_resp_summary/CMAP.up.dn_",pathset.combined,"_",pathset.split,"_",dataset,"_",method,".pdf")
    pdf(file=fname,width=5.5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(4,4,2,2))

  file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_",pathset.combined,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file=file)
  mat.combined <- PATHWAY_CR

  file <- paste0("../output/pathway_conc_resp_summary/old/PATHWAY_CR_",pathset.split,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file=file)
  mat.split <- PATHWAY_CR


  mat.split <- mat.split[is.element(mat.split$pathset,"CMAP"),]
  mat.combined <- mat.combined[is.element(mat.combined$pathset,"CMAP"),]

  mat.split <- mat.split[is.element(mat.split$pathway_class,"estrogen"),]
  mat.combined <- mat.combined[is.element(mat.combined$pathway_class,"estrogen"),]

  dtxsid.list <- c("DTXSID3022536",
                   "DTXSID3037094",
                   "DTXSID5029055",
                   "DTXSID7020182",
                   "DTXSID4022442",
                   "DTXSID8020337",
                   "DTXSID4022369"
  )
  chems <- unique(mat.combined[,c("dtxsid","casrn","name")])
  rownames(chems) <- chems$dtxsid
  chems <- chems[dtxsid.list,]
  nchem <- nrow(chems)
  pathway.list <- unique(mat.combined$pathway)
  npath <- length(pathway.list)

  name.list <- c("dtxsid","casrn","name","pathway","hitcall","bmd","top","hitcall.up","bmd.up","top.up","hitcall.dn","bmd.dn","top.dn")
  res <- as.data.frame(matrix(nrow=npath*nchem,ncol=length(name.list)))
  names(res) <- name.list
  pointer <- 1
  for(dtxsid in dtxsid.list) {
    for(pathway in pathway.list) {
      temp <- mat.combined[is.element(mat.combined$pathway,pathway),]
      temp <- temp[is.element(temp$dtxsid,dtxsid),]
      res[pointer,"dtxsid"] <- temp[1,"dtxsid"]
      res[pointer,"casrn"] <- temp[1,"casrn"]
      res[pointer,"name"] <- temp[1,"name"]
      res[pointer,"pathway"] <- pathway

      res[pointer,"hitcall"] <- temp[1,"hitcall"]
      res[pointer,"bmd"] <- temp[1,"bmd"]
      res[pointer,"top"] <- temp[1,"top"]


      if(res[pointer,"hitcall"] >0.5 && is.na(res[pointer,"top"])) browser()
      pathway.up <- str_replace(pathway,"CMAP","CMAP_up")
      temp <- mat.split[is.element(mat.split$pathway,pathway.up),]
      temp <- temp[is.element(temp$dtxsid,dtxsid),]
      res[pointer,"hitcall.up"] <- temp[1,"hitcall"]
      res[pointer,"bmd.up"] <- temp[1,"bmd"]
      res[pointer,"top.up"] <- temp[1,"top"]

      pathway.dn <- str_replace(pathway,"CMAP","CMAP_dn")
      temp <- mat.split[is.element(mat.split$pathway,pathway.dn),]
      temp <- temp[is.element(temp$dtxsid,dtxsid),]
      res[pointer,"hitcall.dn"] <- temp[1,"hitcall"]
      res[pointer,"bmd.dn"] <- temp[1,"bmd"]
      res[pointer,"top.dn"] <- temp[1,"top"]
      pointer <- pointer+1
    }
  }

  plot(res$top.up~res$top,xlab="top.combined",ylab="top.split",main="CMAP combined vs. split",cex.lab=1.2,cex.axis=1.2,type="n",
       xlim=c(-0.3,0.3),ylim=c(-0.3,0.3))
  lines(c(-10,10),c(0,0))
  lines(c(0,0),c(-10,10))
  lines(c(-10,10),c(-10,10))
  for(i in 1:nrow(res)) {
    x <- res[i,"top"]
    y <- res[i,"top.up"]
    xh <- res[i,"hitcall"]
    yh <- res[i,"hitcall.up"]
    col <- "white"
    if(xh>0.5 && yh>0.5) col= "green"
    if(xh>0.5 && yh<0.5) col= "yellow"
    if(xh<0.5 && yh>0.5) col= "orange"
    points(x,y,pch=21,bg=col)

    x <- res[i,"top"]
    y <- res[i,"top.dn"]
    xh <- res[i,"hitcall"]
    yh <- res[i,"hitcall.dn"]
    col <- "white"
    if(xh>0.5 && yh>0.5) col= "green"
    if(xh>0.5 && yh<0.5) col= "yellow"
    if(xh<0.5 && yh>0.5) col= "orange"
    points(x,y,pch=23,bg=col)
  }
  points(-0.3,0.3,pch=21,bg="green")
  points(-0.3,0.28,pch=21,bg="orange")
  points(-0.3,0.26,pch=21,bg="yellow")
  points(-0.3,0.24,pch=21,bg="white")
  points(-0.3,0.22,pch=21,bg="white")
  points(-0.3,0.20,pch=23,bg="white")

  text(-0.3,0.3,"Hit in both",pos=4,cex=0.9)
  text(-0.3,0.28,"Hit in split",pos=4,cex=0.9)
  text(-0.3,0.26,"Hit in combined",pos=4,cex=0.9)
  text(-0.3,0.24,"Hit in neither",pos=4,cex=0.9)
  text(-0.3,0.22,"up",pos=4,cex=0.9)
  text(-0.3,0.20,"dn",pos=4,cex=0.9)

  if(!to.file)   browser()
  else dev.off()
  file <- paste0("../output/pathway_conc_resp_summary/CMAP.up.dn_",pathset.combined,"_",pathset.split,"_",dataset,"_",method,".xlsx")
  write.xlsx(res,file)
}

