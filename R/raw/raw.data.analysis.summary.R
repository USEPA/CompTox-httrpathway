#source("https://bioconductor.org/biocLite.R")
#.libPaths( "C:/Program Files/R/R-3.3.3/library")
library(fastmatch)
library(reshape2)
#source("R/general_utils.R")
#options(java.parameters = "-Xmx8000m")
#library(httRws)
library(openxlsx)
library(tidyr)
library(grDevices)
library(RColorBrewer)
library(gplots)
library(ape)
#--------------------------------------------------------------------------------------
#' Prepare the raw data for the replicated chemicals
#'
#' @param dir The directory to work in
#'
#--------------------------------------------------------------------------------------
raw.data.analysis.summary <- function(to.file=F,
                                      dir="analysis/raw_count_analysis/") {
  printCurrentFunction()
  if(to.file) {
    welltype <- "dmso"
    fname <- paste0(dir,"raw_data_analysis_summary_heatmap.pdf")
    pdf(file=fname,width=7,height=6,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }

  file <- paste0(dir,"raw_data_analysis_dmso_correlation_mad.xlsx")
  dmso <- read.xlsx(file)
  plate.list <- sort(unique(dmso[,"plate_id"]))
  type.list <- c("bl_dmso","bl_tsa","hbrr","uhrr","TSA","GEN","SIRO","untreated")
  res <- matrix(nrow=length(plate.list),ncol=(length(type.list)+3))
  rownames(res) <- plate.list
  colnames(res) <- c("dmso",type.list,"uhrr_hbrr","bl_dmso_bl_tsa")
  res[] <- 0
  for(plate in plate.list) {
    val <- median(dmso[is.element(dmso[,"plate_id"],plate),"mad"])
    res[plate,"dmso"] <- val
  }

  for(type in type.list) {
    file <- paste0(dir,"raw_data_analysis_l2fc_correlation_mad_",type,".xlsx")
    mat <- read.xlsx(file)
    for(plate in plate.list) {
      val <- median(mat[is.element(mat[,"plate_id"],plate),"mad"])
      res[plate,type] <- val
    }
  }
  file <- paste0(dir,"raw_data_analysis_l2fc_correlation_mad_hbrr_uhrr.xlsx")
  mat <- read.xlsx(file)
  for(plate in plate.list) {
    val <- mat[is.element(mat[,"plate_id"],plate),"mad"]
    res[plate,"uhrr_hbrr"] <- val
  }

  file <- paste0(dir,"raw_data_analysis_l2fc_correlation_mad_bl_dmso_bl_tsa.xlsx")
  mat <- read.xlsx(file)
  for(plate in plate.list) {
    val <- mat[is.element(mat[,"plate_id"],plate),"mad"]
    res[plate,"bl_dmso_bl_tsa"] <- val
  }

  x <- unique(dmso[,c("plate_id","block_id")])
  rownames(x) <- x[,"plate_id"]
  x <- x[rownames(res),]
  col.list <- x[,1]
  col.list[is.element(x[,2],1)] <- "red"
  col.list[is.element(x[,2],2)] <- "black"
  col.list[is.element(x[,2],3)] <- "white"
  col.list[is.element(x[,2],5)] <- "cyan"
  heatmap.2(t(res),
            Rowv=T,
            Colv=T,
            hclustfun=function(x) hclust(d=dist(x),method="ward.D"),
            dendrogram="both",
            symm=F,
            scale="none",
            col=brewer.pal(9,"Reds"),
            margins=c(5,8),
            cexCol=0.1,
            cexRow=1,
            key=T,
            main="MAD values by well type",
            xlab="",
            ylab="",
            trace="none",
            key.title="Key",
            key.xlab="MAD",
            cex.main=1,
            ColSideColors=col.list)
  if(!to.file) browser()
  if(to.file) dev.off()
  file <- paste0(dir,"raw_data_analysis_summary.xlsx")
  mat <- as.data.frame(res)
  mat <- cbind(rownames(mat),mat)
  names(mat)[1] <- "plate_id"
  mat$block_id <- NA
  mat$pg_id <- NA
  name.list <- c("plate_id","pg_id","block_id","dmso","bl_dmso","bl_tsa","hbrr","uhrr","TSA","GEN","SIRO","untreated","uhrr_hbrr","bl_dmso_bl_tsa")
  mat <- mat[,name.list]
  name.list <- c("dmso","bl_dmso","bl_tsa","hbrr","uhrr","TSA","GEN","SIRO","untreated","uhrr_hbrr","bl_dmso_bl_tsa")
  mat$plate_id <- as.character(mat$plate_id)
  x <- mat[,name.list]
  mat$average_mad <- apply(x,FUN=mean,MARGIN=1)
  for(i in 1:nrow(mat)) {
    plate_id <- mat[i,"plate_id"]
    temp <- dmso[is.element(dmso[,"plate_id"],plate_id),]
    block_id <- unique(temp$block_id)
    pg_id <- unique(temp$pg_id)
    mat[i,"pg_id"] <- pg_id
    mat[i,"block_id"] <- block_id
  }

  write.xlsx(mat,file)

}
