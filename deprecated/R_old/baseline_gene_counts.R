#--------------------------------------------------------------------------------------
#' baseline_gene_counts Gene the baseline gene counts for the cell atlas project
#'
#' @param db The name of the Mongo database
#' @param dir The directory where the data will be stored
#'
#' This functions takes files created by export_mongo_httr_well()
#' * httr_cell_atlas
#' * httr_tox21_cpp2
#' @export baseline_gene_counts
#' @importFrom openxlsx write.xlsx
#' @importFrom matrixStats rowMaxs
#' @return nothing, writes out to xlsx
#-------------------------------------------------------------------------------------
baseline_gene_counts <- function(db="httr_cell_atlas",dir="../input/rawdata/cellatlas/") {
  printCurrentFunction() 

  
  file = paste0(dir,db,"_info.RDS")
  info <- readRDS(file)
  file = paste0(dir,db,"_counts.RDS")
  counts <- readRDS(file)

  celltypes = unique(info$trt_name)
  name.list = c("probe","gene",celltypes)
  nprobe = nrow(counts)
  res = as.data.frame(matrix(nrow=nprobe,ncol=length(name.list)))
  names(res) = name.list
  res$probe = rownames(counts)
  probes = rownames(counts)
  z = strsplit(probes,"_")
  df = data.frame(matrix(unlist(z),nrow=length(z),byrow=T))
  res$gene = df[,1]
  for(trt in info$trt_name) {
    sids = info[is.element(info$trt_name,trt),"sample_id"]
    temp = counts[,sids]
    res[,trt] = rowMaxs(temp)
  }
  file = paste0(dir,db,"_baseline_counts.xlsx")
  write.xlsx(res,file)
  for(celltype in celltypes) {
    x = res[,celltype]
    x = x * 1000000 / sum(x)
    res[,celltype] = x
  }
  file = paste0(dir,db,"_baseline_cpm.xlsx")
  write.xlsx(res,file)
}
