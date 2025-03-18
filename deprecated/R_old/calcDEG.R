#--------------------------------------------------------------------------------------
#' calcDEG
#' Calculate the relative variability of genes to get the DEGs
#'
#' @param dataset The name to give to the data set
#' @param dir The directory from which to read all of the raw filesatalog file
#' @param do.read If TRUE, read in the HTTr data file
#' @importFrom openxlsx write.xlsx
#' @importFrom stats mad median
#' @export calcDEG
#' @return nothing
#'
#--------------------------------------------------------------------------------------
calcDEG <- function(dataset="mcf7_ph1_pe1_normal_good_pg",
                    dir="../input/fcdata/",
                    do.read=T) {
  printCurrentFunction(dataset)
  if(do.read) {
    file <- paste0(dir,"CHEM_DICT_",dataset,".RDS")
    print(file)
    CHEM_DICT <- readRDS(file)
    rownames(CHEM_DICT) <- CHEM_DICT$sample_key
    CHEM_DICT <<- CHEM_DICT

    file <- paste0(dir,"FCMAT2_",dataset,".RDS")
    print(file)
    FCMAT2 <- readRDS(file)
    FCMAT2 <- FCMAT2[rownames(CHEM_DICT),]
    FCMAT2 <<- FCMAT2

    cat("  data loaded\n")
  }
  mat <- FCMAT2
  mat[is.na(mat)] <- 0
  chems <- CHEM_DICT
  sk.list <- chems[chems$conc_index<=2,"sample_key"]
  temp <- mat[sk.list,]
  meds <- apply(temp,FUN=median,MARGIN=2)
  mads <- apply(temp,FUN=mad,MARGIN=2)

  temp1 <- sweep(mat,2,meds,"-")
  temp2 <- sweep(temp1,2,mads,"/")

  for(i in c(3,4,5,6)) {
    temp3 <- temp2
    temp3[temp3<i] <- 0
    temp3[temp3>0] <- 1
    rs <- rowSums(temp3,na.rm=T)
    chems <- cbind(chems,rs)
    names(chems)[ncol(chems)] <- paste0("MAD_pos",i)

    temp3 <- -temp2
    temp3[temp3<i] <- 0
    temp3[temp3>0] <- 1
    rs <- rowSums(temp3,na.rm=T)
    chems <- cbind(chems,rs)
    names(chems)[ncol(chems)] <- paste0("MAD_neg",i)

  }

  file <- paste0("../output/gene_deg/",dataset,"_gene_deg.xlsx")
  write.xlsx(chems,file=file)
  file <- paste0("../output/gene_deg/",dataset,"_mad_matrix.RDS")
  saveRDS(temp2, file)
}
