#--------------------------------------------------------------------------------------
#'
#' Calculate PODs for the plier method
#' @param do.laod If TRUE, load the input data into memory
#' @param sigset Name of signature set.
#' @param dataset Name of data set.
#' @param method Pathway scoring method in c("fc", "gsva", "gsea")
#' @param bmr_scale	bmr scaling factor. Default = 1.349
#' @param hccut Remove rows with hitcall less than this value
#'
#' @export
#--------------------------------------------------------------------------------------
plierPOD <- function(do.load=F,
                         sigset="screen_large",
                         dataset="MCF7_pilot_DMEM_6hr_pilot_normal_pe_1",
                         method="gsea",
                         hccut=0.9,
                         cutoff=3,
                         condition="all") {
  printCurrentFunction()
  dir = "../input/httr_mcf7_pilot/fcdata all conditions/"
  dataset.list = c(
    "MCF7_pilot_DMEM_6hr_pilot_normal_pe_1",
    "MCF7_pilot_DMEM_12hr_pilot_normal_pe_1",
    "MCF7_pilot_DMEM_24hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_6hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_12hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_24hr_pilot_normal_pe_1"
  )
  media.list = c("DMEM","DMEM","DMEM","PRF","PRF","PRF")
  time.list = c(6,12,24,6,12,24)

  chems = NULL
  for(dataset in dataset.list) {
    file = paste0(dir,"CHEM_DICT_",dataset,".RData")
    print(file)
    load(file=file)
    chems = rbind(chems,CHEM_DICT)
  }

  dir = "../output/mcf7_pilot/"
  file = paste0(dir,"danilo/TCPL_winning_LVs.xlsx")
  mat = read.xlsx(file)
  mat[is.element(mat$Media,"PRF.DMEM"),"Media"] = "PRF"
  mat[,"name"] = mat[,"Compound_name"]

  mat$dtxsid = NA
  mat$casrn = NA
  for(sid in unique(mat$Compound_ID)) {
    dtxsid = unique(chems[is.element(chems$sample_id,sid),"dtxsid"])
    mat[is.element(mat$Compound_ID,sid),"dtxsid"] = dtxsid
    casrn = unique(chems[is.element(chems$sample_id,sid),"casrn"])
    mat[is.element(mat$Compound_ID,sid),"casrn"] = casrn
  }
  result = unique(mat[,c("dtxsid","casrn","name")])
  result = result[order(result$name),]
  result$dataset = NA
  result$condition = "all"
  result$media = NA
  result$time = NA
  result$nlv = NA
  result$plier_lv_pod_min = NA
  result$plier_lv_pod_min.lci = NA
  result$plier_lv_pod_min.uci = NA
  result$plier_lv_pod_abs5 = NA
  result$plier_lv_pod_abs5.lci = NA
  result$plier_lv_pod_abs5.uci = NA
  result$plier_lv_pod_95 = NA
  result$plier_lv_pod_95.lci = NA
  result$plier_lv_pod_95.uci = NA
  result$plier_lv_burst_pod = NA
  rownames(result) = result$dtxsid

  ###############################################################################################
  # latent variable
  ###############################################################################################
  resall = NULL
  for(i in 1:length(dataset.list)) {
    dataset = dataset.list[i]
    media = media.list[i]
    time = time.list[i]
    result$media = media
    result$time = time
    result$condition = paste(media,time)

    mat1 = mat[is.element(mat$Media,media),]
    mat1 = mat1[mat1$Timepoint==time,]
    for(dtxsid in result$dtxsid) {
      temp = mat1[is.element(mat1$dtxsid,dtxsid),]
      temp = temp[temp$hitcall>hccut,]
      nlv = nrow(temp)
      result[dtxsid,"nlv"] = nlv
      temp = temp[order(temp$bmd),]
      bmdl = temp$bmdl
      bmdu = temp$bmdu
      ratio = bmdu/bmdl
      ratio [is.na(ratio)] = 100
      result[dtxsid,"plier_lv_pod_min"] = temp[1,"bmd"]
      result[dtxsid,"plier_lv_pod_min.lci"] = temp[1,"bmdl"]
      result[dtxsid,"plier_lv_pod_min.uci"] = temp[1,"bmdu"]

      if(nlv>1) {
        bmds = sort(temp$bmd)
        if(nlv>2) {
          db = density(bmds)
          result[dtxsid,"plier_lv_burst_pod"] = db$x[which.max(db$y)]
        }
        minval = 5
        if(minval>nrow(temp)) minval = nrow(temp)
        result[dtxsid,"plier_lv_pod_abs5"] = temp[minval,"bmd"]
        result[dtxsid,"plier_lv_pod_abs5.lci"] = temp[minval,"bmdl"]
        result[dtxsid,"plier_lv_pod_abs5.uci"] = temp[minval,"bmdu"]

        minval = 0.05*nlv
        if(minval<5) minval=5
        if(minval>nrow(temp)) minval = nrow(temp)
        result[dtxsid,"plier_lv_pod_95"] = temp[minval,"bmd"]
        result[dtxsid,"plier_lv_pod_95.lci"] = temp[minval,"bmdl"]
        result[dtxsid,"plier_lv_pod_95.uci"] = temp[minval,"bmdu"]
      }
    }
    result$dataset = dataset
    resall = rbind(resall,result)
  }
  resall[is.na(resall$plier_lv_pod_min),"plier_lv_pod_min"] = 1000
  resall[is.na(resall$plier_lv_pod_95),"plier_lv_pod_95"] = 1000
  resall[is.na(resall$plier_lv_pod_abs5),"plier_lv_pod_abs5"] = 1000
  resall[is.na(resall$plier_lv_burst_pod),"plier_lv_burst_pod"] = 1000

  file = paste0(dir,"plier_lv_pod.xlsx")
  write.xlsx(resall,file)
}

