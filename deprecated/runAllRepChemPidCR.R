#' Run All Replicate Chemical PID Concentration Response
#' 
#' Runs probe ID concentration response for replicate chemicals.
#' 
#' This function has hard-coded dataset names for the replicates. For each
#' replicate, it runs concentration-response directly on the probe ID's. The
#' result is written to disk.
#'
#' @param pval P-value to use for noise estimation. Noise is estimated using two 
#'   lowest concentrations.
#' @param mc.cores Number of cores to use for CR.
#' @param conthits conthits = T uses continuous hitcalls. Continuous hitcalls are
#'   a prerequisitie for using repChemPidPlot(). 
#'
#' @return No output.
#' @export
runAllRepChemPidCR = function(pval = .05, mc.cores = 39, conthits = T){

  #hard-coded dataset names
  studys = c("ph1_", "pilot_")
  floors = c("5","10")
  pes = c("0","1")
  methods = c("none", "normal", "normalold", "apeglm", "ashr")
  combos = expand.grid(studys, floors,pes,methods, stringsAsFactors = F)
  datanames = apply(combos,1,function(x){paste0(c(x,"_pid"),collapse = "")})

  for(dataset in datanames){
    print(dataset)
    #run gene CR
    geneConcResp(dataset=dataset, mc.cores=mc.cores, to.file=T, pval = pval, nametag = NULL, conthits = conthits)
    
  }

}