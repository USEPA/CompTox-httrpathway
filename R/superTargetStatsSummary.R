#--------------------------------------------------------------------------------------
#' Compile the summary statistics for the super targets
#' select chemcial / super targets pairs with good support
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#'
#--------------------------------------------------------------------------------------
superTargetStatsSummary <- function(do.load=F,
                                    sigset="screen_large",
                                    method="fc",
                                    hccut=0.95,
                                    tccut=1.5) {
  printCurrentFunction()
  file = "../input/signatures/signatureDB_master_catalog 2021-03-05.xlsx"
  catalog = read.xlsx(file)
  if(do.load) {
    mat = NULL
    celltype.list = c("MCF7","HepaRG","U2OS")
    dataset.list =c("mcf7_ph1_pe1_normal_block_123","heparg2d_toxcast_pfas_pe1_normal","u2os_toxcast_pfas_pe1_normal")
    for(i in 1:length(celltype.list)) {
      celltype = celltype.list[i]
      dataset = dataset.list[i]
      file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_all.RData")
      print(file)
      load(file=file)
      res.all$celltype = celltype
      mat = rbind(mat,res.all)
    }
    MAT <<- mat
    cmat = unique(mat[,c("sample_id","dtxsid","name","celltype")])
    cmat$burst_bmd = 1000
    cmat$stress_bmd = 1000
    CMAT <<- cmat
  }
  mat = MAT
  mat = mat[mat$bmd_median<1000,]
  cmat = CMAT
  rownames(cmat) = paste(cmat$sample_id,cmat$celltype)
  celltypes = c("MCF7","U2OS","HepaRG")
  sids = unique(cmat$sample_id)
  for(celltype in celltypes) {
    t1 = mat[mat$celltype==celltype,]
    cat(celltype,":",nrow(t1),"\n")
    for(sid in sids) {
      t2 = t1[is.element(t1$sample_id,sid),]
      #if(sid=="EPAPLT0027K06") browser()
      if(nrow(t2)>0) {
        mval = median(t2$bmd_median)
        rn = paste(sid,celltype)
        cmat[rn,"burst_bmd"] = mval
        if(is.element("Stress",t2$super_target)) {
          cmat[rn,"stress_bmd"] = t2[is.element(t2$super_target,"Stress"),"bmd_median"]
        }
      }
    }
  }
  file = paste0("../output/super_target_boxplot/superTargetStatsSummary_CMAT.xlsx")
  write.xlsx(cmat,file)

#browser()
  exclude.list = c("Random","Transcription Factor","Cell Cycle","Biomolecule Process","Cardiovascular")
  mat = mat[!is.element(mat$super_target,exclude.list),]

  set1 = mat[mat$count>=50,]
  mat = mat[mat$count<50,]
  set2 = mat[mat$count>=20,]
  mat = mat[mat$count<20,]
  set3 = mat[mat$count>=10,]
  cat(nrow(set1),nrow(set2),nrow(set3),"\n")
  set1$evidence = "GE.50"
  set2$evidence = "GE.20"
  set3$evidence = "GE.10"
  res = rbind(set1,set2,set3)
  res$burst_bmd = 1000
  res$stress_bmd = 1000
  for(i in 1:nrow(res)) {
    celltype = res[i,"celltype"]
    sid = res[i,"sample_id"]
    rn = paste(sid,celltype)
    res[i,"burst_bmd"] = cmat[rn,"burst_bmd"]
    res[i,"stress_bmd"] = cmat[rn,"stress_bmd"]
  }
  file = paste0("../output/super_target_boxplot/superTargetStatsSummary.xlsx")
  write.xlsx(res,file)
}

