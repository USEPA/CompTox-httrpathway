#----------------------------------------------------------------------------------------
#' Find genes that are often seen with large fold changes
#'
#' @param to.file If TRUE, write the plots to a pdf file
#' @param basedir The base directory to find the raw fold change and chemcial data
#' @param dataset The name of the datset to be used
#'
#' @return nothing is returned
#' @export
#----------------------------------------------------------------------------------------
largeFCgenes<- function(to.file=F,basedir="../input/fcdata/",dataset="DMEM_6hr_pilot_normal_00",l2fc.limit=1.2) {

  file <- paste0(basedir,"FCMAT2_",dataset,".RData")
  print(file)
  load(file)
  file <- paste0(basedir,"CHEM_DICT_",dataset,".RData",sep="")
  load(file)
  rownames(CHEM_DICT) <- CHEM_DICT[,"sample_key"]

  fcmat <- FCMAT2
  chems <- CHEM_DICT
  fcmat[is.na(fcmat)] <- 0

  chems <- chems[order(chems$name,chems$conc),]
  key.list <- chems$sample_key
  fcmat <- fcmat[key.list,]
  result <- NULL
  for(l2fc.limit in c(0.5,0.8,1.0,1.2,1.5,2)) {

    temp <- abs(fcmat)
    temp[temp<l2fc.limit] <- 0
    temp[temp>0] <- 1
    cs <- colSums(temp)
    result <- rbind(result,cs)
    rownames(result)[nrow(result)] <- paste0("l2fc.limit.",l2fc.limit)
  }
  cs <- colSums(result)
  result <- rbind(result,cs)
  rownames(result)[nrow(result)] <- "total"
  result <- t(result)
  result <- result[order(result[,"total"],decreasing=T),]
  file <- paste0("../output/signature_conc_resp_summary/largeFCgenes_ ",dataset," ",l2fc.limit,".xlsx")
  write.xlsx(result,file,row.names=T)
  result <- result[1:50,]

  result <- result[sort(rownames(result)),]
  if(to.file) {
    fname <- paste0("../output/signature_conc_resp_summary/largeFCgenes_ ",dataset," ",l2fc.limit,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  result <- heatmap.2(as.matrix(result),
                      margins=c(10,10),
                      scale="none",
                      main=paste0(dataset),
                      cex.main=0.9,
                      xlab="",
                      ylab="",
                      cexCol=1,
                      cexRow=0.5,
                      col=brewer.pal(8,"Reds"),
                      Rowv=F,
                      Colv=T,
                      trace="none",
                      hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
                      key=T,
                      key.title="Key",
                      key.xlab="counts")

  if(to.file) dev.off()
  else browser()
}
