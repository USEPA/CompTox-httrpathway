#--------------------------------------------------------------------------------------
#' Generate chemicalwise boxplot of the BMD distributions by super_target
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
specificTargetCatalog <- function(method="fc",
                                  celltype="HepaRG",
                                  sigset="screen_large",
                                  hccut=0.95,
                                  tccut=1.5) {
  printCurrentFunction()


  dataset.list = c("heparg2d_toxcast_pfas_pe1_normal","mcf7_ph1_pe1_normal_block_123","u2os_toxcast_pfas_pe1_normal")
  celltype.list = c("HepaRG","MCF7","U2OS")

  res= NULL
  for(i in 1:3) {
    dataset = dataset.list[i]
    celltype = celltype.list[i]
    file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,".xlsx")
    temp = read.xlsx(file)
    res = rbind(res,temp)
  }

  targets = res$specific_targets
  targets = targets[!is.na(targets)]
  targets= targets[!is.element(targets,"-")]

  tlist = NULL
  for(i in 1:length(targets)) {
    t1 = targets[i]
    t2 = str_split(t1,"\\|")[[1]]
    tlist = c(tlist,t2)
  }

  t3 = as.data.frame(table(tlist))
  names(t3) = c("target","frequency")
  file <- paste0("../output/super_target_boxplot/specificTargetCatalog.xlsx")
  write.xlsx(t3,file)
}

