#--------------------------------------------------------------------------------------
#'
#' Summarize the POD overlap with ToxCast
#' @param sigset THe name of the signature set
#' @param dataset Name of the HTTr data set
#' @param method THe signature scoring method
#--------------------------------------------------------------------------------------
signaturePODsummary <- function(sigset="pilot_large_all_100CMAP",
                         dataset="DMEM_6hr_pilot_normal_pe_1",
                         method="gsea") {
  printCurrentFunction()

  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[order(chems$name),]
  dtxsid.list <- chems$dtxsid

  file <- "../toxcast/toxcast_pod.xlsx"
  print(file)
  toxcast <- read.xlsx(file)
  rownames(toxcast) <- toxcast$dtxsid

  file <- "../input/BMDExpress/BMDExpress_Pseudo1_ANOVA_Pathway_Results_Pilot_6h_DMEM.RData"
  load(file=file)
  all_pathway_bmds <- all_pathway_bmds[is.element(all_pathway_bmds$type,"Real"),]

  hc.list = c(0.5,0.6,0.7,0.8,0.9,0.95,0.99)

  name.list = c("hccut","tceq","tclt","tcgt","tceq.count","tclt.count","tcgt.count")
  res = as.data.frame(matrix(nrow=length(hc.list),ncol=length(name.list)))
  names(res) = name.list

  for(i in 1:length(hc.list)) {
    hccut = hc.list[i]
    tceq = 0
    tclt = 0
    tcgt = 0

    tceq.count = 0
    tclt.count = 0
    tcgt.count = 0

    res[i,"hccut"] = hccut
    file <- paste0("../output/signature_pod/signature_pod_",sigset,"_",dataset,"_",method,"_",hccut,".xlsx")

    pod.signature <- read.xlsx(file)
    pod.signature <- pod.signature[order(pod.signature$signature_pod_95),]
    rownames(pod.signature) <- pod.signature$dtxsid
    dtxsid.list <- pod.signature$dtxsid
    toxcast = toxcast[dtxsid.list,]


    for(dtxsid in dtxsid.list) {
      name <- chems[is.element(chems$dtxsid,dtxsid),"name"]

      signature_pod_95 <- pod.signature[dtxsid,"signature_pod_95"]
      signature_pod_95.lci <- pod.signature[dtxsid,"signature_pod_95.lci"]
      signature_pod_95.uci <- pod.signature[dtxsid,"signature_pod_95.uci"]

      signature_pod_95.count <- pod.signature[dtxsid,"signature_pod_95.count"]
      signature_pod_95.count.lci <- pod.signature[dtxsid,"signature_pod_95.count.lci"]
      signature_pod_95.count.uci <- pod.signature[dtxsid,"signature_pod_95.count.uci"]

      toxcast_pod_05 <- toxcast[dtxsid,"pod_uM"]

      #bmds_pod <- min(all_pathway_bmds[is.element(all_pathway_bmds$chem_name,name),"BMD"],na.rm=T)

      if(is.na(signature_pod_95)) signature_pod_95 <- 1000
      if(is.na(signature_pod_95.lci)) signature_pod_95.lci <- 0.001
      if(is.na(signature_pod_95.uci)) signature_pod_95.uci <- 1000

      if(is.na(signature_pod_95.count)) signature_pod_95.count <- 1000
      if(is.na(signature_pod_95.count.lci)) signature_pod_95.count.lci <- 0.001
      if(is.na(signature_pod_95.count.uci)) signature_pod_95.count.uci <- 1000

      if(is.na(toxcast_pod_05)) toxcast_pod_05 <- 1000

      if(signature_pod_95.uci==1000 && signature_pod_95.lci==0.001) signature_pod_95.lci <- 1000
      if(signature_pod_95.count.uci==1000 && signature_pod_95.count.lci==0.001) signature_pod_95.count.lci <- 1000

      delta <- 0.5
      xm <- log10(signature_pod_95.lci)
      x0 <- log10(signature_pod_95)
      xp <- log10(signature_pod_95.uci)
      if(x0-xm < delta) xm <- x0-delta
      if(xp-x0 < delta) xp <- x0+delta
      signature_pod_95.lci <- 10**xm
      signature_pod_95.uci <- 10**xp

      xm <- log10(signature_pod_95.count.lci)
      x0 <- log10(signature_pod_95.count)
      xp <- log10(signature_pod_95.count.uci)
      if(x0-xm < delta) xm <- x0-delta
      if(xp-x0 < delta) xp <- x0+delta
      signature_pod_95.count.lci <- 10**xm
      signature_pod_95.count.uci <- 10**xp

      if(toxcast_pod_05<signature_pod_95.lci) tclt = tclt + 1
      else if(toxcast_pod_05>signature_pod_95.uci) tcgt = tcgt + 1
      else tceq = tceq + 1

      if(toxcast_pod_05<signature_pod_95.count.lci) tclt.count = tclt.count + 1
      else if(toxcast_pod_05>signature_pod_95.count.uci) tcgt.count = tcgt.count + 1
      else tceq.count = tceq.count + 1

    }
    res[i,"tceq"] = tceq
    res[i,"tclt"] = tclt
    res[i,"tcgt"] = tcgt

    res[i,"tceq.count"] = tceq.count
    res[i,"tclt.count"] = tclt.count
    res[i,"tcgt.count"] = tcgt.count
  }
   file <- paste0("../output/signature_pod/signature_pod_summary.xlsx")
  write.xlsx(res,file)
  browser()

}

