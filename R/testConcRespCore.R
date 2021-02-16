testConcRespCore <- function() {
  printCurrentFunction()
  file = "../input/txcplfit2_test.RData"
  load(file=file)

  signaturescoremat = signaturescoremat[1:200]
  for(i in 1:length(signaturescoremat)) {
    x = signaturescoremat[[i]]
    signaturescoremat[[i]]$dtxsid = paste0("dtxsid_",i)
    signaturescoremat[[i]]$sample_id = paste0("sample_id_",i)
    signaturescoremat[[i]]$casrn = paste0("casrn",i)
    signaturescoremat[[i]]$name = paste0("name",i)
  }
  res = NULL
  for(i in 1:length(signaturescoremat)) {
    x = concRespCore(signaturescoremat[[i]],fitmodels,bmd_low_bnd=0.1,bmd_up_bnd=10)
    res = rbind(res,x)
  }

  file = "../input/tcplfit2_test_with_output.RData"
  save(list=c("signaturescoremat","fitmodels","res"),file=file)

  browser()

}
