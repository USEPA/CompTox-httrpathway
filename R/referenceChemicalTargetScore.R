#--------------------------------------------------------------------------------------
#' Generate a quality metric for reference chemicals and targets
#' A chemical that hits its purported target in multiple cell types
#' gets a higher quality score
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#'
#--------------------------------------------------------------------------------------
referenceChemicalTargetScore <- function(do.load=F,
                                         sigset="screen_large",
                                         method="fc",
                                         hccut=0.95,
                                         tccut=1.5) {
  printCurrentFunction()
  if(do.load) {
    res = NULL

    dataset = "heparg2d_toxcast_pfas_pe1_normal"
    celltype = "HepaRG"
    file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_all.RData")
    print(file)
    load(file=file)
    res.all$celltype = celltype
    res = rbind(res,res.all)

    dataset = "mcf7_ph1_pe1_normal_block_123"
    celltype = "MCF7"
    file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_all.RData")
    print(file)
    load(file=file)
    res.all$celltype = celltype
    res = rbind(res,res.all)

    dataset = "u2os_toxcast_pfas_pe1_normal"
    celltype = "U2OS"
    file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_all.RData")
    print(file)
    load(file=file)
    res.all$celltype = celltype
    res = rbind(res,res.all)

    RES <<- res
  }
  mat = RES[RES$match_chem==1,]
  chems = unique(mat[,c("dtxsid","name","use_class","chem_super_target","super_target")])
  chems$MCF7 = -1
  chems$U2OS = -1
  chems$HepaRG = -1
  chems$celltypes = 0
  chems$score = 0
  for(i in 1:nrow(chems)) {
    dtxsid = chems[i,"dtxsid"]
    st = chems[i,"super_target"]
    temp = mat[is.element(mat$dtxsid,dtxsid),]
    temp = temp[is.element(temp$super_target,st),]
    celltypes = 0
    score = 0
    for(j in 1:nrow(temp)) {
      celltype = temp[j,"celltype"]
      active = temp[j,"active"]
      chems[i,celltype] = active
      celltypes = celltypes + 1
      score = score+active
    }
    chems[i,"celltypes"] = celltypes
    chems[i,"score"] = score / (max(1,celltypes))
  }

  file = paste0("../output/super_target_boxplot/referenceChemicalTargetScore.xlsx")
  write.xlsx(chems,file)
}

