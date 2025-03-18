#--------------------------------------------------------------------------------------
#' Normalize FCMAT2 so tht within each plate group, the lowest 2
#' concentrations of each gene have a mean of zero
#'
#' @param chems The CHEMS data frame with chemical information
#' @param method =gene or probe_id
#' @return Global variables are created for the FC matrix (FCMAT2), the SE matrix (SEMAT2)
#' and the chemical dictionary (CHEM_DICT) which translates form the sample key
#' (sample_id_conc_time) to the individual components
#'
#--------------------------------------------------------------------------------------
normalizeFCMAT2 <- function(dataset="DMEM_6hr_screen_normal_pe_1",
                            dir="../input/fcdata/",
                            method="gene",do.read=F)
                             {
  printCurrentFunction(dataset)
  cat("load FCMAT1\n")
  flush.console()
  if(do.read) {
    file <- paste0(dir,"FCMAT2_",dataset,".RData")
    print(file)
    load(file)
    FCMAT2 <<- FCMAT2

    file <- paste0(dir,"CHEM_DICT_",dataset,".RData")
    print(file)
    load(file)
    CHEM_DICT <<- CHEM_DICT

     cat("data loaded\n")
  }
  cat("copy FCMAT2 to mat\n")
  flush.console()
  mat <- FCMAT2
  chems <- CHEM_DICT
  file <- paste0("../input/chemicals/",dataset,"_chemical_map.xlsx")
  chem.map <- read.xlsx(file)
  rownames(chem.map) <- chem.map$sample_key
  res <- NULL
  pg.list <- sort(unique(chem.map$pg_id))
  for(pg in pg.list) {
    sample.list <- chem.map[is.element(chem.map$pg_id,pg),"sample_key"]
    temp <- mat[sample.list,]
    conc.list <- chem.map[sample.list,"conc"]
    mask <- conc.list
    mask[] <- 0
    mask[conc.list<0.05] <- 1
    temp2 <- temp[mask==1,]
    cm <- colMeans(temp2)
    temp3 <- sweep(temp,2,cm,FUN="-")
    res <- rbind(res,temp3)
  }
  file <- paste0(dir,"FCMAT2_",dataset,"_pgnorm.RData")
  FCMAT2 <- res
  save(FCMAT2,file=file)
  file <- paste0(dir,"CHEM_DICT_",dataset,"_pgnorm.RData")
  save(CHEM_DICT,file=file)

}
