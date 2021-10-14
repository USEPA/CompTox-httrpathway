#' Carry out analyses of different POD methods for the pilot study
#'
#'
#' * MCF7_pilot_DMEM_6hr_pilot_normal_pe_1
#' * MCF7_pilot_DMEM_12hr_pilot_normal_pe_1
#' * MCF7_pilot_DMEM_24hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_6hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_12hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_24hr_pilot_normal_pe_1

MCF7PilotPodAnalysis <- function(method="gsea",
                                 celltype="MCF7",
                                 sigset="screen_large",
                                 hccut=0.9,
                                 tccut=1,
                                 cutoff=3,
                                 pval=0.05) {
  printCurrentFunction()
  dir = "../output/mcf7_pilot/"
  nset = 6
  dataset.list = c(
    "MCF7_pilot_DMEM_6hr_pilot_normal_pe_1",
    "MCF7_pilot_DMEM_12hr_pilot_normal_pe_1",
    "MCF7_pilot_DMEM_24hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_6hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_12hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_24hr_pilot_normal_pe_1"
  )
  media.list = c("DMEM","DMEM","DMEM","PRF","PRF","PRF")
  time.list = c(6,12,24,6,12,24)

  # Read in signature results
  mat = NULL
  for(i in 1:nset) {
    dataset = dataset.list[i]
    file = paste0("../output/super_target_boxplot/",celltype,"/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,"_",hccut,"_",tccut,"_",cutoff,"_summary.xlsx")
    temp = read.xlsx(file)
    name.list = c("dtxsid","name","pod_st","pod_sig","nst","nsig")
    temp = temp[,name.list]
    temp1 = temp[,c("dtxsid","name","pod_st","nst")]
    temp1$pod_type = "super_target"
    temp2 = temp[,c("dtxsid","name","pod_sig","nsig")]
    temp2$pod_type = "signature"
    names(temp1)[3] = "pod"
    names(temp2)[3] = "pod"
    names(temp1)[4] = "n"
    names(temp2)[4] = "n"
    temp3 = rbind(temp1,temp2)
    temp3$media = media.list[i]
    temp3$time = time.list[i]
    mat = rbind(mat,temp3)
  }

  for(i in 1:nset) {
    dataset = dataset.list[i]
    file = paste0("../output/gene_conc_resp_summary/GENE_CR_",dataset,"_", pval,"_conthits.RData")
    load(file=file)
    gene = GENE_CR
    dlist = unique(gene$dtxsid)
    name.list = names(mat)
    temp = as.data.frame(matrix(nrow=length(dlist),ncol=length(name.list)))
    names(temp) = name.list
    temp$dtxsid = dlist
    temp$media = media.list[i]
    temp$time = time.list[i]
    temp$pod_type = "gene"
    temp$n = 0
    for(dtxsid in dlist) {
      res = gene[is.element(gene$dtxsid,dtxsid),]
      temp[is.element(temp$dtxsid,dtxsid),"name"] = res[1,"name"]
      res = res[res$hitcall>=hccut,]
      res = res[res$top_over_cutoff>=tccut,]
      blist = res$bmd
      pod = 1000
      if(length(blist)>=10) {
        qlist = quantile(blist,seq(0,1,0.05))
        pod = qlist[2]
      }
      temp[is.element(temp$dtxsid,dtxsid),"pod"] = pod
      temp[is.element(temp$dtxsid,dtxsid),"n"] = length(blist)
    }
    mat = rbind(mat,temp)

  }
  browser()
  file = paste0(dir,"mcf7_pilot_pod_",method,"_",hccut,"_",tccut,"_",cutoff,".xlsx")
  write.xlsx(mat,file)
}
