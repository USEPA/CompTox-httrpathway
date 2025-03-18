#--------------------------------------------------------------------------------------
#' Build a version of the signature concentration-response data with one row
#' per chemical, sample, concentration and one column per signature.
#' The values are top, but set to zero if |top|<cutoff
#'
#'   heparg2d_toxcast_pfas_pe1_normal
#'   mcf7_ph1_pe1_normal_block_123_allPG
#'   u2os_toxcast_pfas_pe1_normal
#'   PFAS_HepaRG
#'   PFAS_U2OS
#'   u2os_pilot_pe1_normal_null_pilot_lowconc
#'   u2os_toxcast_pfas_pe1_normal_refchems
#'   heparg2d_toxcast_pfas_pe1_normal_refchems
#'
#--------------------------------------------------------------------------------------
signatureSliceBuilder <- function(do.load=T,
                                  mc.cores=2,
                                  dataset="mcf7_ph1_pe1_normal_block_123_allPG",
                                  sigset="screen_large",
                                  method="fc",
                                  celltype="MCF7",
                                  hccut=0.95,
                                  tccut=1.5) {
  printCurrentFunction(paste(dataset,sigset,method))


  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    cat(nrow(mat),"\n")
    mat = mat[mat$top_over_cutoff>tccut,]
    cat(nrow(mat),"\n")
    mat = mat[mat$hitcall>hccut,]
    cat(nrow(mat),"\n")
    MAT <<- mat
  }
  mat = MAT
  name.list = c("sample_id","dtxsid","name","signature","cutoff","conc","resp")
  mat = mat[,name.list]

  if(mc.cores > 1){
    cl = makePSOCKcluster(mc.cores)
    result0 = parApply(cl=cl, X=mat, FUN=top.filter, MARGIN = 1,chunk.size = ceiling(nrow(mat)/5/mc.cores) )
    stopCluster(cl)
  }
  else {
    result0 = apply(mat,MARGIN=1,FUN=top.filter)
  }
  result1 = do.call(rbind.data.frame,result0)
  result1$celltype = celltype
  SLICE = result1
  file = paste0("../output/signature_conc_resp_summary/SLICE_",sigset,"_",dataset,"_",method,"_",hccut,"_",tccut,".RData")
  save(SLICE,file=file)
}
top.filter <- function(row) {
  row = as.list(row)
  list2env(row, envir = environment())
  cutoff = as.numeric(cutoff)
  conc = as.numeric(strsplit(conc,"\\|")[[1]])
  resp = as.numeric(strsplit(resp,"\\|")[[1]])
  resp[abs(resp)<cutoff] = 0
  conc = conc[resp!=0]
  resp = resp[resp!=0]
  resp[resp>0] = 1
  resp[resp<0] = -1
  name.list = c("sample_id","dtxsid","name","signature","conc","resp")
  res = as.data.frame(matrix(nrow=length(conc),ncol=length(name.list)))
  names(res) = name.list
  res$sample_id = sample_id
  res$dtxsid = dtxsid
  res$name = name
  res$signature = signature
  res$conc = conc
  res$resp = resp
  return(res)
}

