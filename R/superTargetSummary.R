#--------------------------------------------------------------------------------------
#' Summarize trends in super_target activity
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' mcf7_ph1_pe1_normal_all_pg
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#'
#--------------------------------------------------------------------------------------
superTargetSummary <- function(do.load=F,
                               verbose=F,
                               dataset="u2os_toxcast_pfas_pe1_normal",
                               sigset="screen_large",
                               method="fc",
                               celltype="U2OS",
                               hccut=0.95,
                               tccut=2.5) {
  printCurrentFunction(paste(dataset,sigset,method))

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }
  mat = MAT
  dtxsid.list = unique(MAT$dtxsid)
  mat = mat[mat$hitcall>=hccut,]
  mat = mat[mat$top_over_cutoff>=tccut,]

  st.list = sort(unique(mat$super_target))
  name.list = c("celltype","super_target","nhit","nchem_with_st","nchem_hit","nchem_hit_with_st","nchem_hit_with_st_spec","hit_ratio","spec_ratio")
  res = as.data.frame(matrix(nrow=length(st.list),ncol=length(name.list)))
  names(res) = name.list
  res$celltype = celltype
  file = "../input/chemicals/httr_chemical_annotations 2020-09-14.xlsx"
  rcdb = read.xlsx(file)
  rcdb = rcdb[is.element(rcdb$dtxsid,dtxsid.list),]

  file <- paste0("../output/super_target_boxplot/super_target_boxplot_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  print(file)
  chems = read.xlsx(file)

  for(i in 1:length(st.list)) {
    st = st.list[i]
    cat(st,"\n")
    temp = mat[is.element(mat$super_target,st),]

    ctemp = rcdb[is.element(rcdb$annotation,st),]
    nchem_with_st = length(unique(ctemp$dtxsid))
    nchem_hit_with_st = 0
    nchem_hit_with_st_spec = 0
    ratio1 = 0
    ratio2 = 0
    if(nchem_with_st>0) {
      dlist = ctemp$dtxsid
      temp2 = temp[is.element(temp$dtxsid,dlist),]
      nchem_hit_with_st = length(unique(temp2$dtxsid))
      if(nchem_hit_with_st>0) {
        dlist = unique(temp2$dtxsid)
        count = 0
        for(dtxsid in dlist) {
          if(is.element(dtxsid,chems$dtxsid)) {
            #cat(">>>",dtxsid,"\n")
            temp3 = temp2[is.element(temp2$dtxsid,dtxsid),]
            medval = median(temp3$bmd)
            #print(chems[is.element(chems$dtxsid,dtxsid),])
            burst_pod = min(chems[is.element(chems$dtxsid,dtxsid),"burst_pod"])
            if(medval < burst_pod / 10) count = count+1
          }
        }
        nchem_hit_with_st_spec = count
      }
    }
    if(nchem_with_st) {
      ratio1 = nchem_hit_with_st / nchem_with_st
      ratio2 = nchem_hit_with_st_spec / nchem_with_st
    }
    res[i,"super_target"] = st
    res[i,"nhit"] = nrow(temp)
    res[i,"nchem_with_st"] = nchem_with_st
    res[i,"nchem_hit"] = length(unique(temp$dtxsid))
    res[i,"nchem_hit_with_st"] = nchem_hit_with_st
    res[i,"nchem_hit_with_st_spec"] = nchem_hit_with_st_spec
    res[i,"hit_ratio"] = ratio1
    res[i,"spec_ratio"] = ratio2

    if(verbose && nchem_hit_with_st>0) {
      print(res[i,])
      browser()
    }
  }

  file <- paste0("../output/super_target_boxplot/super_target_summary_",celltype,"_",dataset,"_",sigset,"_",method,".xlsx")
  write.xlsx(res,file)
}

