#--------------------------------------------------------------------------------------
#'
#' Compare the pathway PODs between BMDExpress and DESeq2 / TCPL
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
pathwayBmdeDeseq <- function(to.file=F,
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
    fname <- paste0("../output/pod_laneplot/pathwayBmdeDeseq_",dataset,"_",pathset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))


  file <- "../input/BMDExpress/BMDExpress_Pathway_Results_Pilot_6h_DMEM.RData"
  load(file=file)
  ### all_pathway_bmds

  file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file)
  deseq.all <- PATHWAY_CR
  dtxsid.list <- unique(deseq.all$dtxsid)

  xmin <- 1e-3
  xmax <- 1e3
  res <- NULL
  for(dtxsid in dtxsid.list) {
    name <- chems[is.element(chems$dtxsid,dtxsid),"name"]

    bmds <- all_pathway_bmds[is.element(all_pathway_bmds$chem_name,name),]
    deseq <- deseq.all[is.element(deseq.all$dtxsid,dtxsid),]

      path.list <- sort(unique(bmds$Pathway))
    path.list <- path.list[is.element(path.list,deseq$pathway)]
    rownames(bmds) <- bmds$Pathway
    rownames(deseq) <- deseq$pathway

    bmds <- bmds[path.list,]
    deseq <- deseq[path.list,]

    mask <- deseq$bmd
    mask[] <- 1
    mask[is.na(deseq$bmd)] <- 0

    bmds <- bmds[mask==1,]
    deseq <- deseq[mask==1,]

    b.bmds <- bmds$BMD
    b.deseq <- deseq$bmd

    plot(b.bmds~b.deseq,main=name,cex.axis=1.2,cex.lab=1.2,
         xlab="log(BMD10 DESeq2+tcpl)",ylab="log(BMD10 BMD Express)",
         xlim=c(xmin,xmax),ylim=c(xmin,xmax),log="xy")
    lines(c(1e-5,1e5),c(1e-5,1e5))

    temp <- deseq[,c("dtxsid","casrn","name","pathway","pathway_class","bmdl","bmd","bmdu")]
    temp <- cbind(temp,bmds[,c("chem_name","Pathway","BMDL","BMD","BMDU")])
    res <- rbind(res,temp)
    if(!to.file) browser()
  }
  file <- paste0("../output/pod_laneplot/pathwayBmdeDeseq_",dataset,"_",pathset,"_",method,".xlsx")
  write.xlsx(res,file)
  if(to.file) dev.off()
}

