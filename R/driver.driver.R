#--------------------------------------------------------------------------------------
#' Code to run all calculations
#'
#--------------------------------------------------------------------------------------
driver.driver <- function() {
  printCurrentFunction()

  ####driver(dataset="DMEM_6hr_pilot_normal_00",do.build.random=T)
  ####driver(dataset="DMEM_12hr_pilot_normal_00",do.build.random=T)
  ####driver(dataset="DMEM_24hr_pilot_normal_00",do.build.random=T)
  ####driver(dataset="DMEM_6hr_pilot_none_00",do.build.random=T)
  #driver(dataset="DMEM_12hr_pilot_none_00",do.build.random=T)
  #driver(dataset="DMEM_24hr_pilot_none_00",do.build.random=T)

  driver(dataset="DMEM_6hr_pilot_normal_00",do.run.random=T)
  driver(dataset="DMEM_6hr_pilot_normal_00",do.run.all=T)
  driver(dataset="DMEM_6hr_pilot_normal_00",do.accumulation.plot=T)

  #driver(dataset="DMEM_12hr_pilot_normal_00",do.run.random=T)
  #driver(dataset="DMEM_12hr_pilot_normal_00",do.run.all=T)
  #driver(dataset="DMEM_12hr_pilot_normal_00",do.accumulation.plot=T)

  #driver(dataset="DMEM_24hr_pilot_normal_00",do.run.random=T)
  #driver(dataset="DMEM_24hr_pilot_normal_00",do.run.all=T)
  #driver(dataset="DMEM_24hr_pilot_normal_00",do.accumulation.plot=T)

  #driver(dataset="DMEM_6hr_pilot_none_00",do.run.random=T)
  #driver(dataset="DMEM_6hr_pilot_none_00",do.run.all=T)
  #driver(dataset="DMEM_6hr_pilot_none_00",do.accumulation.plot=T)

  #driver(dataset="DMEM_12hr_pilot_none_00",do.run.random=T)
  #driver(dataset="DMEM_12hr_pilot_none_00",do.run.all=T)
  #driver(dataset="DMEM_12hr_pilot_none_00",do.accumulation.plot=T)

  #driver(dataset="DMEM_24hr_pilot_none_00",do.run.random=T)
  #driver(dataset="DMEM_24hr_pilot_none_00",do.run.all=T)
  #driver(dataset="DMEM_24hr_pilot_none_00",do.accumulation.plot=T)

}
