#--------------------------------------------------------------------------------------
#' Convert the conc-response data to a z score
#' @param do.load If TRUE, load hte large HTTr data set
#' @param mc.cores NUmber of cores to use in multi-core mode=2,
#' @param dataset Name of the HTTr data set being used
#' @param sigset Name of the signature set used
#' @param method Scoring method used
#' @param celltype name of cell type ebing used
#' @param hccut Exclude signature rows with hitcall less than this value
#' @param tccut Exclude signature rows with top_over_cutoff less than this value
#' @importFrom openxlsx read.xlsx
#' @importFrom stats median
#' @export signatureConcRespToZ
#--------------------------------------------------------------------------------------
signatureConcRespToZ <- function(do.load=T,
                                 mc.cores=2,
                                 dataset="heparg2d_toxcast_pfas_pe1_normal",
                                 sigset="screen_large",
                                 method="fc",
                                 celltype="HepaRG",
                                 hccut=0.95,
                                 tccut=1.5) {
  printCurrentFunction(paste(dataset,sigset,method))


  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RDS")
    print(file)
    SIGNATURE_CR <- readRDS(file)
    mat = SIGNATURE_CR
    cat(nrow(mat),"\n")
    mat = mat[mat$top_over_cutoff>tccut,]
    cat(nrow(mat),"\n")
    mat = mat[mat$hitcall>hccut,]
    cat(nrow(mat),"\n")
    MAT <<- mat
  }

  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_summary.xlsx")
  summary = read.xlsx(file)
  mat = MAT
  mat = mat[!is.na(mat$name),]

  zmat = NULL
  slist = unique(mat$sample_id)
  for(sample_id in slist) {
    temp = mat[is.element(mat$sample_id,sample_id),]
    burst = median(temp$bmd)*2
    temp = temp[temp$bmd<burst,]
    temp$z = -(log10(temp$bmd) - log10(burst))*3
    temp$signedz = temp$z * temp$top / abs(temp$top)
    cat(temp[1,"name"],":",range(temp$signedz),"\n")
    if(is.na(temp[1,"name"])) browser()
    zmat = rbind(zmat,temp)
  }
  file = paste0("../output/signature_zmat/ZMAT_",sigset,"_",dataset,"_",method,"_0.05_conthits.RDS")
  ZMAT = zmat
  saveRDS(ZMAT,file)

}
