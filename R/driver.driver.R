#--------------------------------------------------------------------------------------
#' Code to run all calculations
#'
#--------------------------------------------------------------------------------------
driver.driver <- function() {
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
    #driver(dataset=dataset,do.build.random=T)
    #driver(dataset=dataset,do.run.random=T)
    #driver(dataset=dataset,do.run.all=T)
    #driver(dataset=dataset,do.accumulation.plot=T)
    #driver(dataset=dataset,do.pathway.summary.plot=T)
    driver(dataset=dataset,do.pathway.pod=T)
  }

}
