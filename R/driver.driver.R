#--------------------------------------------------------------------------------------
#' Code to run all calculations
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")

#--------------------------------------------------------------------------------------
driver.driver <- function(pathset="PathwaySet_20191107",
                          method="fc",
                          nrandom.chems=1000,
                          mc.cores=30) {
  printCurrentFunction()


  dataset.list = c(
    "DMEM_6hr_pilot_normal_pe_0",
    "DMEM_12hr_pilot_normal_pe_0",
    "DMEM_24hr_pilot_normal_pe_0",

    "DMEM_6hr_pilot_none_pe_0",
    "DMEM_12hr_pilot_none_pe_0",
    "DMEM_24hr_pilot_none_pe_0",

    "DMEM_6hr_pilot_normal_pe_1",
    "DMEM_12hr_pilot_normal_pe_1",
    "DMEM_24hr_pilot_normal_pe_1",

    "DMEM_6hr_pilot_none_pe_1",
    "DMEM_12hr_pilot_none_pe_1",
    "DMEM_24hr_pilot_none_pe_1"#,

    #"PRF_6hr_pilot_normal_pe_0",
    #"PRF_12hr_pilot_normal_pe_0",
    #"PRF_24hr_pilot_normal_pe_0",

    #"PRF_6hr_pilot_none_pe_0",
    #"PRF_12hr_pilot_none_pe_0",
    #"PRF_24hr_pilot_none_pe_0",

    #"PRF_6hr_pilot_normal_pe_1",
    #"PRF_12hr_pilot_normal_pe_1",
    #"PRF_24hr_pilot_normal_pe_1",

    #"PRF_6hr_pilot_none_pe_1",
    #"PRF_12hr_pilot_none_pe_1",
    #"PRF_24hr_pilot_none_pe_1"
  )

  for(dataset in dataset.list) {
    driver(dataset=dataset,pathset=pathset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.build.random=T)
    driver(dataset=dataset,pathset=pathset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.run.random=T)
    driver(dataset=dataset,pathset=pathset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.run.all=T)
    driver(dataset=dataset,pathset=pathset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.accumulation.plot=T)
    driver(dataset=dataset,pathset=pathset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.pathway.summary.plot=T)
    driver(dataset=dataset,pathset=pathset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.pathway.pod=T)
    driver(dataset=dataset,pathset=pathset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.pathway.pod.laneplot=T)
  }

}
