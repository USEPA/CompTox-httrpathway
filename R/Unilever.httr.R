#--------------------------------------------------------------------------------------
#' pull our HTTr data for teh Unilever chemicals
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
Unilever.httr <- function(do.load=F,
                          dataset="u2os_toxcast_pfas_pe1_normal",
                          sigset="screen_large",
                          method="fc",
                          celltype="U2OS",
                          hccut=0.9,
                          tccut=1.5) {
  printCurrentFunction(paste(dataset,sigset,method))

  file = "../input/Unilever/Chemical selection_shortlist_200.xlsx"
  chems = read.xlsx(file)
  dtxsid.list = chems$DTXSID
  nchem = length(dtxsid.list)
  name.list = c("dtxsid","casrn","name","celltype","hasdata","nsigpos","pod")
  res = as.data.frame(matrix(nrow=nchem,ncol=length(name.list)))
  names(res) = name.list

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }

  for(i in 1:nchem) {
    dtxsid = dtxsid.list[i]
    name = chems[[i,"List_CName"]]
    casrn = chems[[i,"CASRN"]]
    cat(dtxsid,"\n")
    res[i,"dtxsid"] = dtxsid
    res[i,"name"] = name
    res[i,"casrn"] = casrn
    res[i,"hasdata"] = "no"
    res[i,"celltype"] = celltype
    if(is.element(dtxsid,MAT$dtxsid)) {
      temp = MAT[is.element(MAT$dtxsid,dtxsid),]
      res[i,"hasdata"] = "yes"
      res[i,"name"] = temp[1,"name"]
      cat(nrow(temp),"\n")
      temp = temp[temp$top_over_cutoff>tccut,]
      cat(nrow(temp),"\n")
      temp = temp[temp$hitcall>hccut,]
      cat(nrow(temp),"\n")
      nsig = nrow(temp)
      res[i,"nsigpos"] = nsig
      res[i,"pod"] = "1000"
      if(nsig>0) {
        blist = temp$bmd
        qb = quantile(blist,probds=seq(0,1,0.05))
        res[i,"pod"]= signif(qb[2],3)
      }
    }
    #if(dtxsid=="DTXSID0023909") browser()
  }
  file = paste0("../input/Unilever/Unilever HTTr overlap ",celltype,".xlsx")
  write.xlsx(res,file)
}

