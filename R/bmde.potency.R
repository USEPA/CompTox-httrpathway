#--------------------------------------------------------------------------------------
#'
#' Export the most potent BMDExpress pathways
#--------------------------------------------------------------------------------------
bmde.potency <- function(to.file=F) {
  printCurrentFunction()

  
  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[order(chems$name),]
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)
  rownames(chems) <- chems$dtxsid
  file <- "../input/BMDExpress/BMDExpress_Pathway_Results_Pilot_6h_DMEM.RData"
  load(file=file)
  bmds <- all_pathway_bmds
  
  allx <- NULL
  ally <- NULL
  labels <- c("10-100","1-10","0.1-1","<0.1")
  temp <- bmds[bmds$BMD>10,]
  y <- temp$N
  x <- y
  x[] <- "A"
  allx <- c(allx,x)
  ally <- c(ally,y)

  temp <- bmds
  temp <- temp[temp$BMD<10,]  
  temp <- temp[temp$BMD>1,]
  y <- temp$N
  x <- y
  x[] <- "B"
  allx <- c(allx,x)
  ally <- c(ally,y)
  
  temp <- bmds
  temp <- temp[temp$BMD<1,]  
  temp <- temp[temp$BMD>0.1,]
  y <- temp$N
  x <- y
  x[] <- "C"
  allx <- c(allx,x)
  ally <- c(ally,y)

  temp <- bmds
  temp <- temp[temp$BMD<0.1,]  
  y <- temp$N
  x <- y
  x[] <- "D"
  allx <- c(allx,x)
  ally <- c(ally,y)
  
  boxplot(ally~allx,cex.lab=1.2,cex.axis=1.2,names=labels,ylab="Number of genes",xlab="BMD (uM)",log="y")
  
  if(!to.file) browser()
  else dev.off()
  
  
  
  
  bmds <- bmds[order(bmds$BMD),]
  bmds <- bmds[bmds$BMD<1,]
  file <- "../input/BMDExpress/BMDExpress potent hits.xlsx"
  #write.xlsx(bmds,file)
 
}

