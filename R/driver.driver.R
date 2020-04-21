#--------------------------------------------------------------------------------------
#' Code to run all calculations
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")

#--------------------------------------------------------------------------------------
driver.driver <- function(sigcatalog="signatureDB_master_catalog 2020-01-31",
                          sigset="pilot_large",
                          nrandom.chems=1000,
                          method="mygsea",
                          mc.cores=30,
                          step1=F,
                          step2=F,
                          step3=F) {
  printCurrentFunction()

  dataset.list = c(
    #"DMEM_6hr_pilot_normal_pe_1"#,
    #"DMEM_12hr_pilot_normal_pe_1",
    #"DMEM_24hr_pilot_normal_pe_1",

    "DMEM_6hr_pilot_normal_pe_0",
    "DMEM_12hr_pilot_normal_pe_0",
    "DMEM_24hr_pilot_normal_pe_0",

    "DMEM_6hr_pilot_none_pe_0",
    "DMEM_12hr_pilot_none_pe_0",
    "DMEM_24hr_pilot_none_pe_0",

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

  #dataset.list <- dataset.list[1]
  # rebuild the FC1MAt and FC2MAT files
  if(step1) driver(mc.cores=mc.cores,do.build.fcmat1.all=T,do.build.fcmat2.all=T)

   # build the random data files
  if(step2) {
    for(dataset in dataset.list) driver(dataset=dataset,sigcatalog=sigcatalog,sigset=sigset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.build.random=T)
  }
  if(step3) {
    for(dataset in dataset.list) {
      tryCatch({
        driver(dataset=dataset,sigcatalog=sigcatalog,sigset=sigset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.run.random=T)
        driver(dataset=dataset,sigcatalog=sigcatalog,sigset=sigset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.run.all=T)
        driver(dataset=dataset,sigcatalog=sigcatalog,sigset=sigset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.signature.summary.plot=T)
        driver(dataset=dataset,sigcatalog=sigcatalog,sigset=sigset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.signature.pod=T)
        driver(dataset=dataset,sigcatalog=sigcatalog,sigset=sigset,method=method,nrandom.chems=nrandom.chems,mc.cores=mc.cores,do.signature.pod.laneplot=T)
      }, warning = function(w) {
      }, error = function(e) {
        cat("ERROR\n")
        print(e)
        browser()
      })
    }
  }
}
