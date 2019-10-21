#' FCMAT for Replicate Chemicals
#' 
#' Generates fold change matrix for the pilot/phase 1 replicate experiment.
#' 
#' Converts deseq2 output to usable FCMAT2 matrices. Also builds CHEM_DICT files
#' as a subset of pre-existing CHEM_DICT files. Generated files are saved
#' directly to disk. Not intended to be run again, but rather to document how
#' it was done originally.
#'
#' @param study Which replicate? "ph1" or "pilot"
#' @param floor Which flooring for input? 5 or 10
#' @param bygene If bygene is TRUE, output will be chemical/concentration by
#'   gene, to be used for pathway analysis; otherwise, it will output 
#'   chemical/concentration by probe id for probe concentration response
#'   modeling.
#'
#' @import data.table
#'
#' @return No output.
#' @export
FCMATrepchems = function(study = "ph1", floor = 10, bygene = T){
  
  # load a master chem dictionary that includes all phase 1 or pilot
  # samples, as necessary.
  if(study == "ph1") load("input/fcdata/CHEM_DICT_unshrunk.RData")
  if(study == "pilot") load("input/fcdata/CHEM_DICT_Pilot44_6.RData")
  masterdict = CHEM_DICT
  
  file = paste0("input/",study,"_repChems_deseq2Output_mean",floor,".RData")
  load(file)
  infile = get(paste0(study,"_repChems"))
  
  chemnames = strsplit(names(infile), "_")
  methodnames = names(infile[[1]])
  
  #access each of the nested lists to separate by processing method
  for(methodname in methodnames){
    bigout = lapply(1:length(infile), function(j){
      bigdata = infile[[j]]
      dtx = chemnames[[j]][1]
      if(study == "ph1") sid = chemnames[[j]][2] else sid = dtx
      data = bigdata[[methodname]]
      conccolnames = colnames(data)[-which(colnames(data) == "probe.id")]
    
      output = lapply(conccolnames, function(conccolname){
        conc = as.numeric(gsub("\\.conc.*","",conccolname))
        mydata = data
        colnames(mydata)[colnames(mydata) == conccolname] = "l2fc"
        if(bygene){
          mydata$probe.id = sub("_.*","",mydata$probe.id)
        }
        mydata = mydata[,list(l2fc = max(l2fc)), by = probe.id]
        
        if(study == "ph1") skey = paste0(sid, "_", conc) else skey = paste0(sid, "_", conc, "_",6)
        mydata$sample_key = skey
        return(mydata)
        
      })
      flatdata = rbindlist(output)
    })
    #unwind the output list and convert to FCMAT2 matrix
    bigflat = rbindlist(bigout)
    FCMAT2 = as.data.frame(dcast(bigflat, sample_key ~ probe.id, value.var = "l2fc", fill = NA_real_))
    rownames(FCMAT2) = FCMAT2$sample_key
    FCMAT2$sample_key = NULL
    FCMAT2 = as.matrix(FCMAT2)
    outname = gsub(" ", "", methodname)
    outname = sub("FPS","", outname)
    if(bygene) finalid = "gene" else finalid = "pid"
    outname = paste0(study,"_", outname,"_",finalid)
    
    #save to disk
    save(FCMAT2, file = paste0("input/fcdata/FCMAT2_", outname,".RData"))
    if(sum(!rownames(FCMAT2) %in% masterdict$sample_key) > 0) warning("some sample keys don't match master dict")
    CHEM_DICT = masterdict[rownames(FCMAT2),]
    
    save(CHEM_DICT, file = paste0("input/fcdata/CHEM_DICT_", outname,".RData"))
    
    
  }
  
  
}