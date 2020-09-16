#--------------------------------------------------------------------------------------
#' Create a candidate list for teh reference chemicals, to be hand annotated
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_toxcast_pfas_pe1_normal
#--------------------------------------------------------------------------------------
refchemdbChemPrep <- function() {
  printCurrentFunction()
  file = "../input/signatures/refchemdb_chem_filtered_unique_with_new.xlsx"
  rcdb = read.xlsx(file)
  rcdb = rcdb[rcdb$refchemdb==1,]
  rcdb = rcdb[,1:3]
  chems = NULL
  dataset.list = c("heparg2d_toxcast_pfas_pe1_normal",
                   "mcf7_ph1_pe1_normal_good_pg",
                   "u2os_toxcast_pfas_pe1_normal")
  chems = NULL
  for(dataset in dataset.list) {
    file = paste0("../input/fcdata/CHEM_DICT_",dataset,".RData")
    load(file=file)
    cd = unique(CHEM_DICT[,c("dtxsid","name")])
    chems = rbind(chems,cd)

  }
  chems = unique(chems)
  chems$super_target = "-"
  rcdb = rbind(rcdb,chems)
  file = "../input/signatures/refchemdb_chem_complete_annotated.xlsx"
  ann = read.xlsx(file)
  dtxsid.list = unique(ann$dtxsid)
  missing = rcdb[!is.element(rcdb$dtxsid,dtxsid.list),]

  file = "../input/signatures/refchemdb_chem_complete.xlsx"
  write.xlsx(rcdb,file)
 file = "../input/signatures/refchemdb_chem_missing.xlsx"
 write.xlsx(missing,file)
}
