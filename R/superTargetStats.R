#--------------------------------------------------------------------------------------
#' Generate hit statistics by super_target
#'
#' @param do.load If TRUE, Load the large input data file
#' @param dataset Name of the HTTr data set
#' @param sigset Name of the signature set
#' @param method Scoring method
#' @param celltype Name of the cell type
#' @param hccut Exclude rows in the data set with hitcall less than this value
#' @param tccut Exclude rows in the data set with top_over_cutoff less than this value
#' @param cutoff The minimum number of signatures hat have to be active in a super
#' target for the super target to be considered active. Default is 5
#' @export
#--------------------------------------------------------------------------------------
superTargetStats <- function(do.load=F,
                             dataset="heparg2d_toxcast_pfas_pe1_normal_refchems",
                             sigset="pilot_tiny",
                             method="fc",
                             celltype="HepaRG",
                             hccut=0.95,
                             tccut=1.5,
                             cutoff=5) {
  printCurrentFunction(paste(dataset,sigset,method))
  if(do.load) {
    file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_",cutoff,"_all.RData")
    print(file)
    load(file=file)
    RES <<- res.all
  }
  mat = RES
  mat[mat$count<2,"active"] = 0
  mat1 = mat[mat$active==1,]
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_active.xlsx")
  write.xlsx(mat1,file)

  st.list = sort(unique(mat$super_target))
  txt = TxT(1,2,3,4,rowname="x")
  name.list = names(txt$mat)
  res = as.data.frame(matrix(nrow=length(st.list),ncol=length(name.list)))
  names(res) = name.list
  for(i in 1:length(st.list)) {
    st = st.list[i]
    temp = mat[is.element(mat$super_target,st),]
    temp0 = temp[temp$active==0,]
    temp1 = temp[temp$active==1,]
    temp00 = temp0[temp0$match_chem==0,]
    temp10 = temp0[temp0$match_chem==1,]
    temp01 = temp1[temp1$match_chem==0,]
    temp11 = temp1[temp1$match_chem==1,]
    a = nrow(temp11)
    b = nrow(temp01)
    c = nrow(temp10)
    d = nrow(temp00)
    txt = TxT(a,b,c,d,rowname=st)
    res[i,] = txt$mat
    if(txt$sens>0.7 && txt$a>=4) {
      cat(st,"\t",txt$a,"\t",txt$sens,"\n")
    }
  }
  file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_",cutoff,"_stats.xlsx")
  write.xlsx(res,file)
}

