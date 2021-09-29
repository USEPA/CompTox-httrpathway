#' Carry out analyses of different POD methods for the pilot study
#' Do hte permutation analyses
#'
#'
#' * MCF7_pilot_DMEM_6hr_pilot_normal_pe_1
#' * MCF7_pilot_DMEM_12hr_pilot_normal_pe_1
#' * MCF7_pilot_DMEM_24hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_6hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_12hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_24hr_pilot_normal_pe_1

MCF7PilotPodPermAnalysis <- function(do.load=F,
                                     dataset="MCF7_pilot_DMEM_6hr_pilot_normal_pe_1",
                                     method="gsea",
                                     celltype="MCF7",
                                     sigset="screen_large",
                                     hccut=0.9,
                                     tccut=1,
                                     cutoff=3,
                                     pval=0.05) {
  printCurrentFunction()
  dir = "../output/mcf7_pilot/"

  file = "../input/signatures/signatureDB_master_catalog 2021-08-27.xlsx"
  catalog = read.xlsx(file)
  catalog = catalog[catalog[,sigset]==1,]
  sslist = paste0(catalog$source,"_",catalog$subsource)
  catalog$ss = sslist
  sslist = unique(sslist)
  scan.grid = expand.grid(seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),
                          seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1),seq(0,1))
  names(scan.grid) = sslist
  rs = rowSums(scan.grid)
  nset = seq(1,nrow(scan.grid),1)
  nperm = 100
  sset = sample(nset,nperm)
  mask = rs
  mask[] = 0
  mask[rs==1] = 1
  mask[rs>=17] = 1
  mask[sset] = 1
  scan.grid = scan.grid[mask==1,]
  ncond = nrow(scan.grid)

  if(do.load) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    MAT <<- SIGNATURE_CR
  }

  res.all = NULL
  for(i in 1:ncond) {
    mat = MAT
    sslist = NULL
    for(j in 1:ncol(scan.grid)) {
      if(scan.grid[i,j]==1) sslist = c(sslist,names(scan.grid)[j])
    }
    cat(paste(sslist,collapse=" "),"\n")
    sig.list = catalog[is.element(catalog$ss,sslist),"parent"]
    nsigtot = length(sig.list)
    mat = mat[is.element(mat$signature,sig.list),]
    if(nrow(mat)>10) {
      result = unique(mat[,c("dtxsid","casrn","name")])
      result = result[order(result$name),]
      result$dataset = dataset
      result$nsig = NA
      result$nsigtot = nsigtot
      result$signature_pod_min = NA
      result$signature_pod_min.lci = NA
      result$signature_pod_min.uci = NA
      result$signature_pod_abs5 = NA
      result$signature_pod_abs5.lci = NA
      result$signature_pod_abs5.uci = NA
      result$signature_pod_95 = NA
      result$signature_pod_95.lci = NA
      result$signature_pod_95.uci = NA
      result$signature_burst_pod = NA
      rownames(result) = result$dtxsid
      sg = scan.grid[1:nrow(result),]
      for(j in 1:nrow(sg)) sg[j,] = scan.grid[i,]
      result = cbind(result,sg)
      ###############################################################################################
      # signature and super_target
      ###############################################################################################
      for(dtxsid in result$dtxsid) {
        temp = mat[is.element(mat$dtxsid,dtxsid),]
        temp = temp[temp$hitcall>hccut,]
        nsig = nrow(temp)
        result[dtxsid,"nsig"] = nsig
        temp = temp[order(temp$bmd),]
        bmdl = temp$bmdl
        bmdu = temp$bmdu
        ratio = bmdu/bmdl
        ratio [is.na(ratio)] = 100
        #temp = temp[ratio<40,]
        result[dtxsid,"signature_pod_min"] = temp[1,"bmd"]
        result[dtxsid,"signature_pod_min.lci"] = temp[1,"bmdl"]
        result[dtxsid,"signature_pod_min.uci"] = temp[1,"bmdu"]

        if(nsig>1) {
          bmds = sort(temp$bmd)
          if(nsig>2) {
            db = density(bmds)
            result[dtxsid,"signature_burst_pod"] = db$x[which.max(db$y)]
          }
          minval = 5
          if(minval>nrow(temp)) minval = nrow(temp)
          result[dtxsid,"signature_pod_abs5"] = temp[minval,"bmd"]
          result[dtxsid,"signature_pod_abs5.lci"] = temp[minval,"bmdl"]
          result[dtxsid,"signature_pod_abs5.uci"] = temp[minval,"bmdu"]

          minval = 0.05*nsig
          if(minval<5) minval=5
          if(minval>nrow(temp)) minval = nrow(temp)
          result[dtxsid,"signature_pod_95"] = temp[minval,"bmd"]
          result[dtxsid,"signature_pod_95.lci"] = temp[minval,"bmdl"]
          result[dtxsid,"signature_pod_95.uci"] = temp[minval,"bmdu"]
        }
      }

      result[is.na(result$signature_pod_min),"signature_pod_min"] = 1000
      result[is.na(result$signature_pod_95),"signature_pod_95"] = 1000
      result[is.na(result$signature_pod_abs5),"signature_pod_abs5"] = 1000
      result[is.na(result$signature_burst_pod),"signature_burst_pod"] = 1000
      res.all = rbind(res.all,result)
    }
  }
  file = paste0(dir,"/signature_perm_pod_",sigset,"_",dataset,"_",method,"_",hccut,"_",cutoff,".xlsx")
  write.xlsx(res.all,file)
}
