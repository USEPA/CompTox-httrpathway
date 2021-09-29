#' Carry out analyses of different POD methods for the pilot study
#'
#'
#' * MCF7_pilot_DMEM_6hr_pilot_normal_pe_1
#' * MCF7_pilot_DMEM_12hr_pilot_normal_pe_1
#' * MCF7_pilot_DMEM_24hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_6hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_12hr_pilot_normal_pe_1
#' * MCF7_pilot_PRF_24hr_pilot_normal_pe_1

MCF7PilotPodAnalysisDebug <- function(method="gsea",
                                 celltype="MCF7",
                                 sigset="screen_large",
                                 hccut=0.9,
                                 tccut=1,
                                 cutoff=3,
                                 pval=0.05) {
  printCurrentFunction()
  dir = "../output/mcf7_pilot/"
  dataset = "MCF7_pilot_DMEM_6hr_pilot_normal_pe_1"
  file = paste0("../input/fcdata/FCMAT2_",dataset,".RData")
  load(file=file)
  fnew = FCMAT2

  dataset = "DMEM_6hr_pilot_normal_pe_1"
  file = paste0("../input/fcdata/FCMAT2_",dataset,".RData")
  load(file=file)
  fold = FCMAT2

  browser()
  dataset = "MCF7_pilot_DMEM_6hr_pilot_normal_pe_1"
  file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file=file)
  mat.new = SIGNATURE_CR


  #SIGNATURE_CR_pilot_large_all_100CMAP_DMEM_6hr_pilot_normal_pe_1_mygsea_0.05_conthits
  method = "mygsea"
  sigset = "pilot_large_all_100CMAP"
  dataset = "DMEM_6hr_pilot_normal_pe_1"
  file = paste0("../output/signature_conc_resp_summary/pilot/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
  print(file)
  load(file=file)
  mat.old = SIGNATURE_CR

  sigs = unique(mat.old$signature)
  sigs = sigs[is.element(sigs,mat.new$signature)]

  mat.new = mat.new[is.element(mat.new$signature,sigs),]
  mat.old = mat.old[is.element(mat.old$signature,sigs),]

  sig = "CMAP estradiol 1e-07 100 8281 100"
  dtxsid = "DTXSID5029055"
  tnew = mat.new[is.element(mat.new$dtxsid,dtxsid),]
  tnew = tnew[is.element(tnew$signature,sig),]

  told = mat.old[is.element(mat.old$dtxsid,dtxsid),]
  told = told[is.element(told$signature,sig),]
  browser()
}
