#--------------------------------------------------------------------------------------
#'
#' Get the pathway ranks for chemicals
#--------------------------------------------------------------------------------------
pathwayRank <- function(to.file=F,
                        dataset="DMEM_6hr_pilot_normal_pe_1",
                        pathset="PathwaySet_20191107",
                        method="mygsea") {
  printCurrentFunction()

  if(to.file) {
    fname <- paste0("../output/pathway_rank/pathway_rank_",dataset,"_",pathset,"_",method,".pdf")
    pdf(file=fname,width=5.5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(4,4,2,2))

  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[order(chems$name),]
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)
  rownames(chems) <- chems$dtxsid
  file <- "../input/BMDExpress/BMDExpress_Pathway_Results_Pilot_6h_DMEM.RData"
  load(file=file)
  bmds <- all_pathway_bmds

  file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file)
  deseq <- PATHWAY_CR
  drank <- abs(deseq$hitcall * log10(deseq$bmd10))
  deseq$drank <- drank
  deseq <- deseq[order(deseq$drank,decreasing=T),]


  name.list <- c("dtxsid","casrn","name","target_key",
                 "deseq2.first.pathway","deseq2.first.pathway.score","deseq2.first.pathway.class",
                 "deseq.first.ontarget.pathway","deseq.first.ontarget.pathway.rank",
                 "deseq.first.ontarget.pathway.score","deseq.first.ontarget.pathway.pod",
                 "bmde.first.pathway","bmde.first.pathway.class",
                 "bmde.first.ontarget.pathway","bmde.first.ontarget.pathway.rank","bmde.first.ontarget.pathway.pod")
  result <- as.data.frame(matrix(nrow=nchem,ncol=length(name.list)))
  names(result) <- name.list
  result$dtxsid <- chems$dtxsid
  rownames(result) <- result$dtxsid

  for(dtxsid in dtxsid.list) {
    temp <- deseq[is.element(deseq$dtxsid,dtxsid),]
    result[dtxsid,"name"] <- temp[1,"name"]
    result[dtxsid,"casrn"] <- temp[1,"casrn"]
    target_key <- chems[dtxsid,"target_key"]
    result[dtxsid,"target_key"] <- target_key
    cat(result[dtxsid,"name"],":",target_key,"\n")

    pathset.list <- NULL

    if(target_key=="thyroid") pathset.list <- c("thyroid")
    if(target_key=="ER") pathset.list <- c("estrogen","Tamoxifen","Cyproterone")
    if(target_key=="ion channel") pathset.list <- c("ion channel","Amiodarone")
    if(target_key=="PPAR") pathset.list <- c("ppar")
    if(target_key=="CYPs") pathset.list <- c("p450","conazole")
    if(target_key=="AR") pathset.list <- c("Testosterone","androgen","Nilutamide","Flutamide")
    if(target_key=="mitochondria") pathset.list <-c("mitochondria")
    if(target_key=="cholesterol") pathset.list <- c("statin","sterol processing","steroid synthesis","glitazone","fibrate","adipogenesis","fatty acid")

    if(target_key=="electron chain") pathset.list <- NULL
    if(target_key=="protein reactive") pathset.list <- NULL
    if(target_key=="adrenergic") pathset.list <- NULL
    if(target_key=="protein synthesis") pathset.list <- NULL
    if(target_key=="Plant PPO") pathset.list <- NULL
    if(target_key=="DNA") pathset.list <- NULL

    if(!is.null(pathset.list)) {
      #DESEQ
      mask <- vector(length=nrow(temp),mode="integer")
      mask[] <- 0
      mask[is.element(temp$pathway_class,pathset.list)] <- 1
      result[dtxsid,"deseq2.first.pathway"] <- temp[1,"pathway"]
      result[dtxsid,"deseq2.first.pathway.score"] <- temp[1,"drank"]
      result[dtxsid,"deseq2.first.pathway.class"] <- temp[1,"pathway_class"]
      index <- min(which(mask==1))

      result[dtxsid,"deseq.first.ontarget.pathway"] <- temp[index,"pathway"]
      result[dtxsid,"deseq.first.ontarget.pathway.rank"] <- index
      result[dtxsid,"deseq.first.ontarget.pathway.score"] <- temp[index,"drank"]
      result[dtxsid,"deseq.first.ontarget.pathway.pod"] <- temp[index,"bmd10"]

      #BMDS
      temp2 <- bmds[is.element(bmds$chem_name,result[dtxsid,"name"]),]
      rownames(temp2) <- temp2$Pathway
      rownames(temp) <- temp$pathway
      plist <- temp2$Pathway
      temp3 <- temp[plist,"pathway_class"]
      temp2$pathway_class <-  temp3
      temp2 <- temp2[order(temp2$BMDL),]
      mask <- vector(length=nrow(temp2),mode="integer")
      mask[] <- 0
      mask[is.element(temp2$pathway_class,pathset.list)] <- 1
      result[dtxsid,"bmde.first.pathway"] <- temp2[1,"Pathway"]
      result[dtxsid,"bmde.first.pathway.class"] <- temp2[1,"pathway_class"]
      if(sum(mask)>0) {
        index <- min(which(mask==1))

        result[dtxsid,"bmde.first.ontarget.pathway"] <- temp2[index,"Pathway"]
        result[dtxsid,"bmde.first.ontarget.pathway.rank"] <- index
        result[dtxsid,"bmde.first.ontarget.pathway.pod"] <- temp2[index,"BMD"]
      }
    }
  }
  file <-  file <- "../input/pathway_dictionary.xlsx"
  dict <- read.xlsx(file)

  plot(c(1,1),type="n",xlim=c(1,2000),ylim=c(1,2000),cex.axis=1.2,cex.lab=1.2,
       #main="DESeq2 vs BMDExpress Target Pathway Rank",
       main=paste(dataset,":",method),
       xlab="DESeq2 rank",ylab="BMDExpress rank",
       log="xy")
  lines(c(1,10000),c(1,10000))
  for(val in c(5,10,50,100,500,1000)) {
    lines(c(1,2000),c(val,val),col="gray")
    lines(c(val,val),c(1,2000),col="gray")
  }
  for(i in 1:nrow(result)) {
    if(!is.na(result[i,"deseq.first.ontarget.pathway.rank"]) && !is.na(result[i,"bmde.first.ontarget.pathway.rank"])) {
      x <- result[i,"deseq.first.ontarget.pathway.rank"]
      x <- x+rnorm(1,0,x*0.05)
      y <- result[i,"bmde.first.ontarget.pathway.rank"]
      target_key <- result[i,"target_key"]
      if(target_key=="thyroid") pathset.list <- c("thyroid")
      if(target_key=="ER") pathset.list <- c("estrogen","Tamoxifen","Cyproterone")
      if(target_key=="ion channel") pathset.list <- c("ion channel","Amiodarone")
      if(target_key=="PPAR") pathset.list <- c("ppar")
      if(target_key=="CYPs") pathset.list <- c("p450","conazole")
      if(target_key=="AR") pathset.list <- c("Testosterone","androgen","Nilutamide","Flutamide")
      if(target_key=="mitochondria") pathset.list <-c("mitochondria")
      if(target_key=="cholesterol") pathset.list <- c("statin","sterol processing","steroid synthesis","glitazone","fibrate","adipogenesis","fatty acid")

      color <- dict[is.element(dict$pathway_class,pathset.list[1]),"color"]
      points(x,y,bg=color,pch=21,cex=2)
    }
  }

  plot(c(1,1),type="n",xlim=c(1e-4,1e3),ylim=c(1e-4,1e3),cex.axis=1.2,cex.lab=1.2,
       main=paste(dataset,":",method),
       xlab="DESeq2 POD (uM)",ylab="BMDExpress POD (uM)",
       log="xy")
  lines(c(0.000001,10000),c(0.000001,10000))
  for(val in c(0.001,0.01,0.1,1,10,100)) {
    lines(c(1e-5,1e5),c(val,val),col="gray")
    lines(c(val,val),c(1e-5,1e5),col="gray")
  }

  for(i in 1:nrow(result)) {
    if(!is.na(result[i,"deseq.first.ontarget.pathway.pod"]) && !is.na(result[i,"bmde.first.ontarget.pathway.pod"])) {
      x <- result[i,"deseq.first.ontarget.pathway.pod"]
      y <- result[i,"bmde.first.ontarget.pathway.pod"]
      target_key <- result[i,"target_key"]
      if(target_key=="thyroid") pathset.list <- c("thyroid")
      if(target_key=="ER") pathset.list <- c("estrogen","Tamoxifen","Cyproterone")
      if(target_key=="ion channel") pathset.list <- c("ion channel","Amiodarone")
      if(target_key=="PPAR") pathset.list <- c("ppar")
      if(target_key=="CYPs") pathset.list <- c("p450","conazole")
      if(target_key=="AR") pathset.list <- c("Testosterone","androgen","Nilutamide","Flutamide")
      if(target_key=="mitochondria") pathset.list <-c("mitochondria")
      if(target_key=="cholesterol") pathset.list <- c("statin","sterol processing","steroid synthesis","glitazone","fibrate","adipogenesis","fatty acid")

      color <- dict[is.element(dict$pathway_class,pathset.list[1]),"color"]
      points(x,y,bg=color,pch=21,cex=2)
    }
  }
  if(!to.file) browser()
  else dev.off()
  file <- paste0("../output/pathway_rank/pathway_rank_",dataset,"_",pathset,"_",method,".xlsx")
  write.xlsx(result,file)
}

