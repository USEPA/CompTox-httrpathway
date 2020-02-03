#--------------------------------------------------------------------------------------
#'
#' Get the pathway ranks for chemicals
#--------------------------------------------------------------------------------------
pathwayRank.DESEQ2 <- function(to.file=F,
                               dataset1="DMEM_6hr_pilot_none_pe_1",
                               dataset2="DMEM_6hr_pilot_normal_pe_1",
                               pathset="PathwaySet_20191107",
                               method="gsva",
                               cutoff=0.5) {
  printCurrentFunction()

  if(to.file) {
    fname <- paste0("../output/pathway_rank/pathway_rank_",dataset1,"_",dataset2,"_",pathset,"_",method,".pdf")
    pdf(file=fname,width=5.5,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(2,1),mar=c(4,4,2,2))

  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[order(chems$name),]
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)
  rownames(chems) <- chems$dtxsid

  file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset1,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file)
  deseq1 <- PATHWAY_CR
  drank <- abs(deseq1$hitcall * log10(deseq1$bmd))
  deseq1$drank <- drank
  deseq1 <- deseq1[order(deseq1$drank,decreasing=T),]

  file <- paste0("../output/pathway_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset2,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file)
  deseq2 <- PATHWAY_CR
  drank <- abs(deseq2$hitcall * log10(deseq2$bmd))
  deseq2$drank <- drank
  deseq2 <- deseq2[order(deseq2$drank,decreasing=T),]

  name.list <- c("dtxsid","casrn","name","target_key",
                 "deseq1.first.pathway","deseq1.first.pathway.score","deseq1.first.pathway.class",
                 "deseq1.first.ontarget.pathway","deseq1.first.ontarget.pathway.rank",
                 "deseq1.first.ontarget.pathway.score","deseq1.first.ontarget.pathway.pod",
                 "deseq1.first.ontarget.pathway.top_over_cutoff",
                 "deseq1.pathway.hits",

                 "deseq2.first.pathway","deseq2.first.pathway.score","deseq2.first.pathway.class",
                 "deseq2.first.ontarget.pathway","deseq2.first.ontarget.pathway.rank",
                 "deseq2.first.ontarget.pathway.score","deseq2.first.ontarget.pathway.pod",
                 "deseq2.first.ontarget.pathway.top_over_cutoff",
                 "deseq2.pathway.hits"
                 )
  result <- as.data.frame(matrix(nrow=nchem,ncol=length(name.list)))
  names(result) <- name.list
  result$dtxsid <- chems$dtxsid
  rownames(result) <- result$dtxsid

  for(dtxsid in dtxsid.list) {
    #DESEQ 1
    temp <- deseq1[is.element(deseq1$dtxsid,dtxsid),]
    result[dtxsid,"name"] <- temp[1,"name"]
    result[dtxsid,"casrn"] <- temp[1,"casrn"]
    target_key <- chems[dtxsid,"target_key"]
    result[dtxsid,"target_key"] <- target_key
    cat(result[dtxsid,"name"],":",target_key,"\n")
    result[dtxsid,"deseq1.pathway.hits"] <- nrow(temp[temp$hitcall>cutoff,])
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
      mask <- vector(length=nrow(temp),mode="integer")
      mask[] <- 0
      mask[is.element(temp$pathway_class,pathset.list)] <- 1
      result[dtxsid,"deseq1.first.pathway"] <- temp[1,"pathway"]
      result[dtxsid,"deseq1.first.pathway.score"] <- temp[1,"drank"]
      result[dtxsid,"deseq1.first.pathway.class"] <- temp[1,"pathway_class"]
      index <- min(which(mask==1))

      result[dtxsid,"deseq1.first.ontarget.pathway"] <- temp[index,"pathway"]
      result[dtxsid,"deseq1.first.ontarget.pathway.rank"] <- index
      result[dtxsid,"deseq1.first.ontarget.pathway.score"] <- temp[index,"drank"]
      result[dtxsid,"deseq1.first.ontarget.pathway.pod"] <- temp[index,"bmd"]
      result[dtxsid,"deseq1.first.ontarget.pathway.top_over_cutoff"] <- temp[index,"top_over_cutoff"]
    }

    #DESEQ 2
    temp <- deseq2[is.element(deseq2$dtxsid,dtxsid),]
    result[dtxsid,"name"] <- temp[1,"name"]
    result[dtxsid,"casrn"] <- temp[1,"casrn"]
    target_key <- chems[dtxsid,"target_key"]
    result[dtxsid,"target_key"] <- target_key
    cat(result[dtxsid,"name"],":",target_key,"\n")
    result[dtxsid,"deseq2.pathway.hits"] <- nrow(temp[temp$hitcall>cutoff,])

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
      mask <- vector(length=nrow(temp),mode="integer")
      mask[] <- 0
      mask[is.element(temp$pathway_class,pathset.list)] <- 1
      result[dtxsid,"deseq2.first.pathway"] <- temp[1,"pathway"]
      result[dtxsid,"deseq2.first.pathway.score"] <- temp[1,"drank"]
      result[dtxsid,"deseq2.first.pathway.class"] <- temp[1,"pathway_class"]
      index <- min(which(mask==1))

      result[dtxsid,"deseq2.first.ontarget.pathway"] <- temp[index,"pathway"]
      result[dtxsid,"deseq2.first.ontarget.pathway.rank"] <- index
      result[dtxsid,"deseq2.first.ontarget.pathway.score"] <- temp[index,"drank"]
      result[dtxsid,"deseq2.first.ontarget.pathway.pod"] <- temp[index,"bmd"]
      result[dtxsid,"deseq2.first.ontarget.pathway.top_over_cutoff"] <- temp[index,"top_over_cutoff"]

    }
  }
  file <-  file <- "../input/pathway_dictionary.xlsx"
  dict <- read.xlsx(file)
  ##########################################################################################
  # Rank plot
  ##########################################################################################
  plot(c(1,1),type="n",xlim=c(1,2000),ylim=c(1,2000),cex.axis=1.2,cex.lab=1.2,
       main="DESeq2 Normal vs. No Shrinkage Pathway Rank",
       xlab="No Shrinkage",ylab="Normal Shrinkage",
       log="xy")
  lines(c(1,10000),c(1,10000))
  for(val in c(5,10,50,100,500,1000)) {
    lines(c(1,2000),c(val,val),col="gray")
    lines(c(val,val),c(1,2000),col="gray")
  }
  for(i in 1:nrow(result)) {
    if(!is.na(result[i,"deseq1.first.ontarget.pathway.rank"]) && !is.na(result[i,"deseq2.first.ontarget.pathway.rank"])) {
      x <- result[i,"deseq1.first.ontarget.pathway.rank"]
      x <- x+rnorm(1,0,x*0.05)
      y <- result[i,"deseq2.first.ontarget.pathway.rank"]
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

  ##########################################################################################
  # Potency plot
  ##########################################################################################
  plot(c(1,1),type="n",xlim=c(1e-4,1e3),ylim=c(1e-4,1e3),cex.axis=1.2,cex.lab=1.2,
       main="DESeq2 Normal vs. No Shrinkage Pathway POD",
       xlab="No Shrinkage POD (uM)",ylab="Normal Shrinkage POD (uM)",
       log="xy")
  lines(c(0.000001,10000),c(0.000001,10000))
  for(val in c(0.001,0.01,0.1,1,10,100)) {
    lines(c(1e-5,1e5),c(val,val),col="gray")
    lines(c(val,val),c(1e-5,1e5),col="gray")
  }

  for(i in 1:nrow(result)) {
    if(!is.na(result[i,"deseq1.first.ontarget.pathway.pod"]) && !is.na(result[i,"deseq2.first.ontarget.pathway.pod"])) {
      x <- result[i,"deseq1.first.ontarget.pathway.pod"]
      y <- result[i,"deseq2.first.ontarget.pathway.pod"]
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

  ##########################################################################################
  # Top over cutoff plot
  ##########################################################################################
  plot(c(1,1),type="n",xlim=c(0,5),ylim=c(0,5),cex.axis=1.2,cex.lab=1.2,
       main="DESeq2 Normal vs. No Shrinkage Pathway POD",
       xlab="No Shrinkage top/cutoff",ylab="Normal Shrinkage top/cutoff")
  lines(c(0,10000),c(0,10000))
  for(val in c(0,1,2,3,4,5)) {
    lines(c(1e-5,1e5),c(val,val),col="gray")
    lines(c(val,val),c(1e-5,1e5),col="gray")
  }

  for(i in 1:nrow(result)) {
    if(!is.na(result[i,"deseq1.first.ontarget.pathway.pod"]) && !is.na(result[i,"deseq2.first.ontarget.pathway.pod"])) {
      x <- result[i,"deseq1.first.ontarget.pathway.top_over_cutoff"]
      y <- result[i,"deseq2.first.ontarget.pathway.top_over_cutoff"]
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

  ##########################################################################################
  # Pahtway hits plot
  ##########################################################################################
  plot(c(1,1),type="n",xlim=c(0,1000),ylim=c(0,1000),cex.axis=1.2,cex.lab=1.2,
       main="DESeq2 Normal vs. No Shrinkage Pathway POD",
       xlab="No Shrinkage Pathway Hits",ylab="Normal Shrinkage Pathway Hits")
  lines(c(0,10000),c(0,10000))
  for(val in c(0,200,400,600,800,1000)) {
    lines(c(1e-5,1e5),c(val,val),col="gray")
    lines(c(val,val),c(1e-5,1e5),col="gray")
  }

  for(i in 1:nrow(result)) {
    if(!is.na(result[i,"deseq1.first.ontarget.pathway.pod"]) && !is.na(result[i,"deseq2.first.ontarget.pathway.pod"])) {
      x <- result[i,"deseq1.pathway.hits"]
      y <- result[i,"deseq2.pathway.hits"]
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
  #file <- paste0("../output/pathway_rank/pathway_rank_",dataset,"_",pathset,"_",method,".xlsx")
  #write.xlsx(result,file)
}

