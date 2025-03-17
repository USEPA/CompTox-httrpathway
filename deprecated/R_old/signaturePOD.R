#--------------------------------------------------------------------------------------
#'
#' Calculate PODs at the signature level
#' @param do.load If TRUE, load the input data into memory
#' @param sigset Name of signature set.
#' @param dataset Name of data set.
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param hccut Remove rows with hitcall less than this value
#' @importFrom openxlsx write.xlsx
#' @importFrom stats density median
#'
#' @export signaturePOD
#--------------------------------------------------------------------------------------
signaturePOD <- function(do.load=F,
                         sigset="screen_large",
                         dataset="MCF7_pilot_DMEM_6hr_pilot_normal_pe_1",
                         method="gsea",
                         hccut=0.9,
                         cutoff=3,
                         condition="all") {
  printCurrentFunction()

  if(do.load) {
    file <- paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RDS")
    print(file)
    MAT <<- readRDS(file)
    file <- paste0("../output/gene_conc_resp_summary/GENE_CR_",dataset,"_0.05_conthits.RDS")
    print(file)
    GMAT <<- readRDS(file)
  }
  mat = MAT
  gmat = GMAT
  result = unique(mat[,c("dtxsid","casrn","name")])
  result = result[order(result$name),]
  result$dataset = dataset
  result$condition = condition
  result$nsig = NA
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
  result$nst = NA
  result$super_target_pod = NA
  result$super_target_burst_pod = NA
  result$ngene = NA
  result$gene_pod_min = NA
  result$gene_pod_min.lci = NA
  result$gene_pod_min.uci = NA
  result$gene_pod_abs5 = NA
  result$gene_pod_abs5.lci = NA
  result$gene_pod_abs5.uci = NA
  result$gene_pod_95 = NA
  result$gene_pod_95.lci = NA
  result$gene_pod_95.uci = NA
  result$gene_burst_pod = NA
  rownames(result) = result$dtxsid

  ###############################################################################################
  # gene
  ###############################################################################################
  for(dtxsid in result$dtxsid) {
    temp = gmat[is.element(gmat$dtxsid,dtxsid),]
    temp = temp[temp$hitcall>hccut,]
    ngene = nrow(temp)
    result[dtxsid,"ngene"] = ngene
    temp = temp[order(temp$bmd),]
    bmdl = temp$bmdl
    bmdu = temp$bmdu
    ratio = bmdu/bmdl
    ratio [is.na(ratio)] = 100
    temp = temp[ratio<40,]
    result[dtxsid,"gene_pod_min"] = temp[1,"bmd"]
    result[dtxsid,"gene_pod_min.lci"] = temp[1,"bmdl"]
    result[dtxsid,"gene_pod_min.uci"] = temp[1,"bmdu"]

    if(ngene>1) {
      bmds = sort(temp$bmd)
      if(ngene>2) {
        db = density(bmds)
        result[dtxsid,"gene_burst_pod"] = db$x[which.max(db$y)]
      }
      minval = 5
      if(minval>nrow(temp)) minval = nrow(temp)
      result[dtxsid,"gene_pod_abs5"] = temp[minval,"bmd"]
      result[dtxsid,"gene_pod_abs5.lci"] = temp[minval,"bmdl"]
      result[dtxsid,"gene_pod_abs5.uci"] = temp[minval,"bmdu"]

      minval = 0.05*ngene
      if(minval<5) minval=5
      if(minval>nrow(temp)) minval = nrow(temp)
      result[dtxsid,"gene_pod_95"] = temp[minval,"bmd"]
      result[dtxsid,"gene_pod_95.lci"] = temp[minval,"bmdl"]
      result[dtxsid,"gene_pod_95.uci"] = temp[minval,"bmdu"]
    }
  }

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

      ##################################################################################
      # super_target
      ##################################################################################
      st.list = unique(temp$super_target)
      st.list = st.list[!is.na(st.list)]
      if(length(st.list)>0) {
        st.pods = as.data.frame(matrix(nrow=length(st.list),ncol=2))
        names(st.pods) = c("super_target","pod")
        st.pods[,1] = st.list
        st.pods[,2] = NA
        rownames(st.pods) = st.pods[,1]
        for(st in st.list) {
          tempst = temp[is.element(temp$super_target,st),]
          if(nrow(tempst)>=cutoff) st.pods[st,"pod"] = median(tempst$bmd)
        }
        st.pods = st.pods[!is.na(st.pods$pod),]
        if(nrow(st.pods)>0) {
          st.pods = st.pods[order(st.pods$pod),]

          result[dtxsid,"nst"] = nrow(st.pods)
          result[dtxsid,"super_target_pod"] = st.pods[1,"pod"]
          bmds = sort(st.pods$pod)
          if(nrow(st.pods)>2) {
            db = density(bmds)
            result[dtxsid,"super_target_burst_pod"] = db$x[which.max(db$y)]
          }
        }
      }
    }
  }

  result[is.na(result$signature_pod_min),"gene_pod_min"] = 1000
  result[is.na(result$signature_pod_95),"gene_pod_95"] = 1000
  result[is.na(result$signature_pod_abs5),"gene_pod_abs5"] = 1000
  result[is.na(result$signature_burst_pod),"gene_burst_pod"] = 1000

  result[is.na(result$signature_pod_min),"signature_pod_min"] = 1000
  result[is.na(result$signature_pod_95),"signature_pod_95"] = 1000
  result[is.na(result$signature_pod_abs5),"signature_pod_abs5"] = 1000
  result[is.na(result$signature_burst_pod),"signature_burst_pod"] = 1000

  result[is.na(result$super_target_pod),"super_target_pod"] = 1000
  result[is.na(result$super_target_burst_pod),"super_target_burst_pod"] = 1000
  result[is.na(result$nst),"nst"] = 0

  file = paste0("../output/signature_pod/signature_pod_",condition,"_",sigset,"_",dataset,"_",method,"_",hccut,"_",cutoff,".xlsx")
  write.xlsx(result,file)
}

