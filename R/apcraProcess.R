#--------------------------------------------------------------------------------------
#' Build th4e APCRA data sets
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_toxcast_pfas_pe1_normal
#--------------------------------------------------------------------------------------
apcraProcess <- function(do.load=T,
                         sigset="screen_large",
                         method="fc",
                         hccut=0.95) {
  printCurrentFunction()
  file = "../input/chemicals/apcra_ret_pro_list.xlsx"
  temp = read.xlsx(file)
  dtxsid.list = unique(temp$DTXSID)

  file = paste0("../output/signature_refchemdb/validated_signatures_merged.xlsx")
  temp=read.xlsx(file)
  sig.list = temp$signature

  ds.list = c("heparg2d_toxcast_pfas_pe1_normal",
              "mcf7_ph1_pe1_normal_good_pg",
              "u2os_toxcast_pfas_pe1_normal")
  ct.list  = c("HepaRG","MCF7","U2OS")
  for(i in 1:3) {
    dataset = ds.list[i]
    celltype = ct.list[i]
    cat(dataset,":",celltype,"\n")
    if(do.load) {
      file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
      print(file)
      load(file=file)
      MAT <<- SIGNATURE_CR
    }
    mat = MAT
    cat(nrow(mat),"\n")
    mat = mat[is.element(mat$dtxsid,dtxsid.list),]
    cat(nrow(mat),"\n")
    mat = mat[mat$hitcall>=hccut,]
    cat(nrow(mat),"\n")
    SIGNATURE_CR = mat
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_APCRA_",celltype,"_hit_only_",method,"_0.05_conthits.RData")
    save(SIGNATURE_CR,file=file)
    mat = mat[is.element(mat$signature,sig.list),]
    cat(nrow(mat),"\n")
    SIGNATURE_CR = mat
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_APCRA_",celltype,"_hit_only_filtered_",method,"_0.05_conthits.RData")
    save(SIGNATURE_CR,file=file)
  }
}

