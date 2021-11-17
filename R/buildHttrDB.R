#--------------------------------------------------------------------------------------
#' BUild the mysql database for the httr data
#'

#--------------------------------------------------------------------------------------
buildHttrDB <- function(do.clean=F,
                        do.init=F,
                        do.load.datasets=F,
                        do.load.pods=F,
                        sigcatalog="signatureDB_master_catalog 2021-09-29",
                        sigset="screen_large",
                        method="gsea",
                        hccut=0.9,
                        tccut=1) {
  printCurrentFunction(paste(sigset,method))
  db = "res_httr"
  dslist = c(
    "MCF7_pilot_DMEM_6hr_pilot_normal_pe_1",
    "MCF7_pilot_DMEM_12hr_pilot_normal_pe_1",
    "MCF7_pilot_DMEM_24hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_6hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_12hr_pilot_normal_pe_1",
    "MCF7_pilot_PRF_24hr_pilot_normal_pe_1",
    "heparg2d_toxcast_pfas_pe1_normal_v2",
    "mcf7_ph1_pe1_normal_block_123_allPG",
    "u2os_toxcast_pfas_pe1_normal_v2",
    "u2os_pilot_pe1_normal_null_pilot",
    "u2os_toxcast_pfas_pe1_normal_v2_refchems",
    "heparg2d_toxcast_pfas_pe1_normal_v2_refchems",
    "tox21_cpp5_u2os_pe1_normal",
    "tox21_cpp5_heparg_pe1_normal"
  )
  if(do.clean) {
    runQuery("delete from concresp",db)
    runQuery("delete from pod",db)
    runQuery("delete from sample",db)
    runQuery("delete from chemical",db)
    runQuery("delete from signature",db)
    runQuery("delete from signature_list",db)
    runQuery("delete from study",db)
  }
  if(do.init) {
    file = "mysql/dataset_index.xlsx"
    mat = read.xlsx(file)
    mat$plate_effect = as.character(mat$plate_effect)
    runQuery("delete from study",db)
    runInsertTable(mat,"study",db=db,do.halt=T,verbose=F,get.id=T)

    file = paste0("../input/signatures/signatureDB.RData")
    load(file=file)
    rownames(sigdb) = sigdb$signature

    file = paste0("../input/signatures/",sigcatalog,".xlsx")
    catalog = read.xlsx(file)
    rownames(catalog) = catalog$signature

    sigs = sigdb$signature
    sigs = sigs[is.element(sigs,catalog$signature)]
    sigdb = sigdb[sigs,]
    catalog = catalog[sigs,]
    sigdb$super_target = catalog$super_target
    sigdb$effect_direction = catalog$effect_direction
    sigdb$target_class = catalog$target_class
    sigdb$super_target_level = catalog$super_target_level
    nlist = names(sigdb)
    nlist[is.element(nlist,"gene.list")] = "genelist"
    names(sigdb) = nlist
    runQuery("delete from signature",db)
    runInsertTable(sigdb,"signature",db=db,do.halt=T,verbose=F,get.id=T)

    catalog = catalog[!is.na(catalog[,sigset]),]

    temp2 = catalog[,c("parent","signature")]
    temp2[,1] = sigset
    names(temp2)[1] = "sigset"
    runQuery("delete from signature_list",db)
    runInsertTable(temp2,"signature_list",db=db,do.halt=T,verbose=F,get.id=T)
  }
  if(do.load.datasets) {

    nds = length(dslist)
    #nds = 1
    for(i in 1:nds) {
      dataset = dslist[i]
      cat(">>> load conc-resp data for dataset:",dataset,"\n")
      file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
      print(file)
      load(file=file)
      mat = SIGNATURE_CR
      mat = mat[mat$hitcall>hccut,]
      mat = mat[mat$top_over_cutoff>tccut,]
      chems = unique(mat[,c("dtxsid","casrn","name")])
      samps = unique(mat[,c("dtxsid","sample_id")])
      runInsertTable(chems,"chemical",db=db,do.halt=T,verbose=F,get.id=T)
      runInsertTable(samps,"sample",db=db,do.halt=T,verbose=F,get.id=T)
      nchem = runQuery("select count(*) from chemical",db)[1,1]
      nsamp = runQuery("select count(*) from sample",db)[1,1]
      cat("chems, samples:",nchem,nsamp,"\n")
      fields = runQuery("desc concresp",db)[,"Field"]
      nlist = names(mat)
      nlist = nlist[is.element(nlist,fields)]
      temp = mat[,nlist]
      temp$dataset = dataset
      temp$method = method
      temp$sigset = sigset
      runQuery(paste0("delete from concresp where dataset='",dataset,"'"),db)
      cat("size the of dataset:",nrow(temp),"\n")
      runInsertTable(temp,"concresp",db=db,do.halt=T,verbose=F,get.id=T)
    }
  }
  if(do.load.pods) {
    nds = length(dslist)
    nds = 6
    for(i in 1:nds) {
      dataset = dslist[i]
      cat(">>> load pods for dataset:",dataset,"\n")
      condition = "all"
      hccut = 0.9
      cutoff = 3
      file = paste0("../output/signature_pod/signature_pod_",condition,"_",sigset,"_",dataset,"_",method,"_",hccut,"_",cutoff,".xlsx")
      mat = read.xlsx(file)

      tlist = c("signature_pod_min","signature_pod_abs5","signature_pod_95",
                "gene_pod_min","gene_pod_abs5","gene_pod_95")

      runQuery(paste0("delete from pod where dataset='",dataset,"'"),db)
      for(type in tlist) {
        uci = paste0(type,".uci")
        lci = paste0(type,".lci")
        temp = mat[,c("dtxsid",type,uci,lci)]
        names(temp) = c("dtxsid","pod","pod_uci","pod_lci")
        temp$pod_method = type
        temp$dataset = dataset
        temp$method = method
        runInsertTable(temp,"pod",db=db,do.halt=T,verbose=F,get.id=T)
      }
    }
  }

}

