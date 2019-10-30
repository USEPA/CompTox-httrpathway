library(openxlsx)
#--------------------------------------------------------------------------------------
#' Code to run all calculations
#'
#--------------------------------------------------------------------------------------
driver <- function(dataset="DMEM_6hr_pilot_normal_00",
                   pathset="PathwaySet_20191031",
                   nrandom.chems=1000,
                   mc.cores=30,
                   do.build=F,
                   do.build.random=F,
                   do.run.random=F,
                   do.run.all=F,
                   do.accumulation.plot=F,
                   do.all=F) {
  printCurrentFunction()

  if(do.build) {
    buildFCMAT1(dataset=dataset,dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_normal/DMEM_6/",filetype="tsv")
    buildFCMAT2(dataset=dataset,dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
  }
  if(do.build.random) {
    randomdata(dataset=dataset, nchem=nrandom.chems)
  }

  nullset <- paste0(dataset,"_RAND",nrandom.chems)
  if(do.run.random){
    runAllPathwayCR_pval(dataset=nullset,
                         nullset=nullset,
                         pathset=pathset,
                         method = "fc",
                         do.plot = F,
                         mc.cores = c(mc.cores,mc.cores))
  }
  if(do.run.all){
    runAllPathwayCR_pval(dataset=dataset,
                         nullset=nullset,
                         pathset=pathset,
                         method = "fc",
                         do.plot = T,
                         mc.cores = c(mc.cores,mc.cores))
    cat("Look for output in \n
        ../output/pathway_score_summary/\n
        ../output/pathway_conc_resp_plots/\n
        ../output/pathway_conc_resp_summary/\n
        \n")
  }
  if(do.accumulation.plot) {
    pathwayAccumNullPlot(pathset=pathset,
                         dataset=dataset,
                         method="fc",
                         nullset=nullset,
                         nametag="conthits",
                         mc.cores=mc.cores)
    cat("Look for output in ../output/accumplots/ \n")
  }
}
