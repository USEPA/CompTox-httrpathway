#--------------------------------------------------------------------------------------
#'
#' Build lane plots by chemical list and signature class, across the datasets
#' @param to.file If TRUE, write plots to a file
#--------------------------------------------------------------------------------------
signatureClassHM <- function(to.file=F,
                           dataset="DMEM_6hr_pilot_normal_pe_1",
                           pathset="PathwaySet_20191107",
                           method = "gsva",
                           threshold=0.5) {
  printCurrentFunction()
  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)

  if(to.file) {
    fname <- paste0("../output/pod_laneplot/signatureClassHM_",dataset,"_",pathset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  file <- paste0("../output/signature_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file=file)
  mat <- PATHWAY_CR

  pathclass.list <- sort(unique(mat$signature_class))

  pathclass.list <- c(
    "Amiodarone",
    "androgen",
    "apoptosis",
    "cell cycle",
    "conazole",
    "Cyproterone",
    "cytotoxicity",
    "dna damage",
    "estrogen",
    "fibrate",
    "Flutamide",
    "glitazone",
    "gpcr",
    "hdac",
    "hypoxia",
    "ion channel",
    "microtubule",
    "mitochondria",
    "nfkb",
    "Nilutamide",
    "other",
    "oxidative stress",
    "p450",
    "ppar",
    "random",
    "statin",
    "steroid synthesis",
    "sterol processing",
    "Tamoxifen",
    "Testosterone",
    "thyroid",
    "tnf"
  )

  npathclass <- length(pathclass.list)

  for(pathclass in pathclass.list) {
    temp <- mat[is.element(mat$signature_class,pathclass),]
    n <- length(unique(temp$signature))
    cat(pathclass,n,"\n")
  }

  res <- as.data.frame(matrix(nrow=nchem,ncol=npathclass))
  rownames(res) <- dtxsid.list
  names(res) <- pathclass.list
  res[] <- 0
  for(dtxsid in dtxsid.list) {
    temp1 <- mat[is.element(mat$dtxsid,dtxsid),]
    for(pathclass in pathclass.list) {
      temp2 <- temp1[is.element(temp1$signature_class,pathclass),]
      if(nrow(temp2)>0) {
        bot <- nrow(temp2)
        top <- nrow(temp2[temp2$hitcall>threshold,])
        res[dtxsid,pathclass]  <- top/bot
      }
    }
  }
  rownames(res) <- chems$name

  result <- heatmap.2(as.matrix(res),
                      margins=c(10,10),
                      scale="none",
                      main=paste(dataset,"\n",method,":",pathset),
                      xlab="",
                      ylab="",
                      cexCol=0.7,
                      cexRow=0.7,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="l2fc",
                      cex.main=1)

  if(!to.file) browser()

  # low conc version
  res <- as.data.frame(matrix(nrow=nchem,ncol=npathclass))
  rownames(res) <- dtxsid.list
  names(res) <- pathclass.list
  res[] <- 0
  for(dtxsid in dtxsid.list) {
    temp1 <- mat[is.element(mat$dtxsid,dtxsid),]
    for(pathclass in pathclass.list) {
      temp2 <- temp1[is.element(temp1$signature_class,pathclass),]
      if(nrow(temp2)>0) {
        bot <- nrow(temp2)
        mask <- temp2$hitcall
        mask[] <- 1
        mask[temp2$hitcall<threshold] <- 0
        mask[temp2$bmd>10] <- 0
        top <- sum(mask)
        res[dtxsid,pathclass]  <- top/bot
      }
    }
  }
  rownames(res) <- chems$name

  result <- heatmap.2(as.matrix(res),
                      margins=c(10,10),
                      scale="none",
                      main=paste(dataset,"\n",method,":",pathset),
                      xlab="",
                      ylab="",
                      cexCol=0.7,
                      cexRow=0.7,
                      Rowv=T,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      col=brewer.pal(9,"Reds"),
                      key.title="Key",
                      key.xlab="l2fc",
                      cex.main=1)

  if(!to.file) browser()

    if(to.file) dev.off()
}

