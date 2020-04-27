library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
#--------------------------------------------------------------------------------------
#' Code to run all calculations
#' @param method signature scoring method in c("fc", "gsva", "mygsea")
#--------------------------------------------------------------------------------------
driver <- function(dataset="DMEM_6hr_pilot_normal_pe_1",
                   sigcatalog="signatureDB_master_catalog 2020-04-05",
                   sigset="pilot_large_all_100CMAP",
                   nrandom.chems=1000,
                   mc.cores=30,
                   method="mygsea",
                   do.build.fcmat1.all=F,
                   do.build.fcmat2.all=F,
                   do.build.random=F,
                   do.run.random=T,
                   do.run.all=T,
                   do.signature.summary.plot=T,
                   do.signature.pod=T,
                   do.signature.pod.laneplot=T,
                   do.all=F) {
  printCurrentFunction(paste(dataset,":",sigset))


  if(do.build.fcmat1.all) {
    buildFCMAT1(dataset="DMEM_6hr_pilot_normal_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_normal/DMEM_6/",filetype="tsv")
    buildFCMAT1(dataset="DMEM_12hr_pilot_normal_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_normal/DMEM_12/",filetype="tsv")
    buildFCMAT1(dataset="DMEM_24hr_pilot_normal_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_normal/DMEM_24/",filetype="tsv")

    buildFCMAT1(dataset="DMEM_6hr_pilot_normal_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_normal/DMEM_6/",filetype="tsv")
    buildFCMAT1(dataset="DMEM_12hr_pilot_normal_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_normal/DMEM_12/",filetype="tsv")
    buildFCMAT1(dataset="DMEM_24hr_pilot_normal_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_normal/DMEM_24/",filetype="tsv")

    buildFCMAT1(dataset="DMEM_6hr_pilot_none_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_none/DMEM_6/",filetype="tsv")
    buildFCMAT1(dataset="DMEM_12hr_pilot_none_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_none/DMEM_12/",filetype="tsv")
    buildFCMAT1(dataset="DMEM_24hr_pilot_none_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_none/DMEM_24/",filetype="tsv")

    buildFCMAT1(dataset="DMEM_6hr_pilot_none_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_none/DMEM_6/",filetype="tsv")
    buildFCMAT1(dataset="DMEM_12hr_pilot_none_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_none/DMEM_12/",filetype="tsv")
    buildFCMAT1(dataset="DMEM_24hr_pilot_none_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_none/DMEM_24/",filetype="tsv")

    #buildFCMAT1(dataset="PRF_6hr_pilot_normal_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_normal/PRF_6/",filetype="tsv")
    #buildFCMAT1(dataset="PRF_12hr_pilot_normal_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_normal/PRF_12/",filetype="tsv")
    #buildFCMAT1(dataset="PRF_24hr_pilot_normal_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_normal/PRF_24/",filetype="tsv")

    #buildFCMAT1(dataset="PRF_6hr_pilot_normal_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_normal/PRF_6/",filetype="tsv")
    #buildFCMAT1(dataset="PRF_12hr_pilot_normal_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_normal/PRF_12/",filetype="tsv")
    #buildFCMAT1(dataset="PRF_24hr_pilot_normal_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_normal/PRF_24/",filetype="tsv")

    #buildFCMAT1(dataset="PRF_6hr_pilot_none_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_none/PRF_6/",filetype="tsv")
    #buildFCMAT1(dataset="PRF_12hr_pilot_none_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_none/PRF_12/",filetype="tsv")
    #buildFCMAT1(dataset="PRF_24hr_pilot_none_pe_0",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_0-shrinkage_none/PRF_24/",filetype="tsv")

    #buildFCMAT1(dataset="PRF_6hr_pilot_none_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_none/PRF_6/",filetype="tsv")
    #buildFCMAT1(dataset="PRF_12hr_pilot_none_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_none/PRF_12/",filetype="tsv")
    #buildFCMAT1(dataset="PRF_24hr_pilot_none_pe_1",dir="../input/httr_mcf7_pilot/meanncnt0_5-plateteffect_1-shrinkage_none/PRF_24/",filetype="tsv")
  }

  if(do.build.fcmat2.all) {
    buildFCMAT2(dataset="DMEM_6hr_pilot_normal_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    buildFCMAT2(dataset="DMEM_12hr_pilot_normal_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    buildFCMAT2(dataset="DMEM_24hr_pilot_normal_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")

    buildFCMAT2(dataset="DMEM_6hr_pilot_none_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    buildFCMAT2(dataset="DMEM_12hr_pilot_none_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    buildFCMAT2(dataset="DMEM_24hr_pilot_none_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")

    buildFCMAT2(dataset="DMEM_6hr_pilot_normal_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    buildFCMAT2(dataset="DMEM_12hr_pilot_normal_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    buildFCMAT2(dataset="DMEM_24hr_pilot_normal_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")

    buildFCMAT2(dataset="DMEM_6hr_pilot_none_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    buildFCMAT2(dataset="DMEM_12hr_pilot_none_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    buildFCMAT2(dataset="DMEM_24hr_pilot_none_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")

    #buildFCMAT2(dataset="PRF_6hr_pilot_normal_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    #buildFCMAT2(dataset="PRF_12hr_pilot_normal_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    #buildFCMAT2(dataset="PRF_24hr_pilot_normal_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")

    #buildFCMAT2(dataset="PRF_6hr_pilot_none_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    #buildFCMAT2(dataset="PRF_12hr_pilot_none_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    #buildFCMAT2(dataset="PRF_24hr_pilot_none_pe_0",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")

    #buildFCMAT2(dataset="PRF_6hr_pilot_normal_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    #buildFCMAT2(dataset="PRF_12hr_pilot_normal_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    #buildFCMAT2(dataset="PRF_24hr_pilot_normal_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")

    #buildFCMAT2(dataset="PRF_6hr_pilot_none_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    #buildFCMAT2(dataset="PRF_12hr_pilot_none_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
    #buildFCMAT2(dataset="PRF_24hr_pilot_none_pe_1",dir="../input/fcdata/",method="gene",do.read=T,chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx")
  }
  if(do.build.random) {
    randomdata(dataset=dataset, nchem=nrandom.chems)
  }

  nullset <- paste0(dataset,"_RAND",nrandom.chems)
  if(do.run.random || do.all){
    runAllSignatureCR(dataset=nullset,
                      nullset=nullset,
                      sigset=sigset,
                      sigcatalog=sigcatalog,
                      method = method,
                      do.plot = F,
                      mc.cores = c(mc.cores,mc.cores))
  }
  if(do.run.all || do.all){
    runAllSignatureCR(dataset=dataset,
                      nullset=nullset,
                      sigset=sigset,
                      sigcatalog=sigcatalog,
                      method = method,
                      do.plot = T,
                      mc.cores = c(mc.cores,mc.cores))
    cat("Look for output in \n
        ../output/signature_score_summary/\n
        ../output/signature_conc_resp_plots/\n
        ../output/signature_conc_resp_summary/\n
        \n")
  }
  if(do.signature.summary.plot || do.all) {
    signatureClassSummaryPlot(to.file=T,dataset=dataset,
                              sigcatalog=sigcatalog,
                              sigset=sigset,
                              method = method)
  }
  if(do.signature.pod || do.all) {
    signaturePOD(sigset=sigset,
                 dataset=dataset,
                 method=method,
                 hit.threshold=0.5)
  }
  if(do.signature.pod.laneplot) {
    podLaneplot(to.file=T,
                dataset=dataset,
                sigset=sigset,
                method=method)
  }
}
