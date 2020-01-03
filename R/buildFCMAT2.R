#--------------------------------------------------------------------------------------
#' Transpose and filter the fold change matrix FCMAT1 in long format into a
#' gene x sample format.
#'
#' @param chems The CHEMS data frame with chemical information
#' @param method =gene or probe_id
#' @return Global variables are created for the FC matrix (FCMAT2), the SE matrix (SEMAT2)
#' and the chemical dictionary (CHEM_DICT) which translates form the sample key
#' (sample_id_conc_time) to the individual components
#'
#--------------------------------------------------------------------------------------
buildFCMAT2 <- function(dataset="DMEM_6hr_pilot_normal_pe_0",dir="../input/fcdata/",
                        method="gene",do.read=T,
                        chemical.file="../input/chemicals/HTTr.Sample.Matrix.2017.04.24.xlsx") {
  printCurrentFunction()
  if(do.read) {
    file <- paste0(dir,"FCMAT1_",dataset,".RData")
    print(file)
    load(file)
    FCMAT1 <<- FCMAT1
    cat("data loaded\n")
  }
  if(!exists("DSSTOX")) {
    file <- "../input/chemicals/DSSTox_2019-05-16.RData"
    load(file=file)
    DSSTOX <<- DSSTOX
  }
  mat <- FCMAT1
  dsstox <- DSSTOX
  rownames(dsstox) <- dsstox$dsstox_substance_id
  chem.map <- read.xlsx(chemical.file)
  sid.list <- unique(mat$sample_key)
  temp <- str_split(sid.list,"_")
  name.list <- c("sample_key","sample_id","conc","time","casrn","name","dtxsid","media","conc_index")
  chems <- as.data.frame(matrix(nrow=length(sid.list),ncol=length(name.list)))
  names(chems) <- name.list
  chems$sample_key <- sid.list
  for(i in 1:nrow(chems)) {
    row <- temp[[i]]
    conc_index <- row[2]
    sid <- row[1]
    time <- row[5]
    time <- as.numeric(str_replace(time,"h",""))
    conc <- row[3]
    conc <- as.numeric(str_replace(conc,"uM",""))
    temp1 <- chem.map[is.element(chem.map$EPA_Sample_ID,sid),]
    dtxsid <- temp1[1,"DTXSID"]
    casrn <- dsstox[dtxsid,"casrn"]
    name <- temp1[1,"Chem.Name"]
    media <- row[4]
    chems[i,"sample_id"] <- sid
    chems[i,"conc"] <- conc
    chems[i,"time"] <- time
    chems[i,"casrn"] <- casrn
    chems[i,"name"] <- name
    chems[i,"dtxsid"] <- dtxsid
    chems[i,"media"] <- media
    chems[i,"conc_index"] <- conc_index
  }

  CHEM_DICT <- chems
  file <- paste0(dir,"CHEM_DICT_",dataset,".RData")
  save(CHEM_DICT,file=file)

  if(method=="gene") {
	  mat2 <- mat[,c("sample_key","gene","l2fc")]
	  matfc <- acast(mat2,gene~sample_key,mean)

	  matp <- mat[,c("sample_key","gene","padj")]
	  matpv <- acast(matp,gene~sample_key,mean)
	  matpv <- t(matpv)

	  mats <- mat[,c("sample_key","gene","se")]
	  matse <- acast(mats,gene~sample_key,mean)
	  matse <- t(matse)
  }
	else if(method=="probe_id") {
	  mat2 <- mat[,c("sample_key","probe_id","l2fc")]
	  matfc <- acast(mat2,probe_id~sample_key,mean)

	  matp <- mat[,c("sample_key","probe_id","padj")]
	  matpv <- acast(matp,probe_id~sample_key,mean)
	  matpv <- t(matpv)

	  mats <- mat[,c("sample_key","probe_id","se")]
	  matse <- acast(mats,probe_id~sample_key,mean)
	  matse <- t(matse)

	}


  mat3 <- matfc
	mat3[is.nan(mat3)] <- NA
	cat("Dimension of mat3 before filtering: ",dim(mat3),"\n")
	temp <- mat3
	temp[!is.na(mat3)] <- 1
	temp[is.na(mat3)] <- 0
	rs <- rowSums(temp)
	nsample <- dim(mat3)[2]
	mask <- rs
	mask[] <- 1
	mask[rs<0.95*nsample] <- 0
	matfc <- matfc[mask==1,]
	cat("Dimension of mat3 after filtering: ",dim(matfc),"\n")
	matfc <- t(matfc)

  FCMAT2 <- matfc
  file <- paste0(dir,"FCMAT2_",dataset,".RData")
  save(FCMAT2,file=file)
  print(dim(FCMAT2))

  matpv <- matpv[,colnames(FCMAT2)]
  matse <- matse[,colnames(FCMAT2)]
  matse.inv <- 1/matse
  matfcse <- matse.inv * FCMAT2

  FCMAT2.FCoverSE <- matfcse
  print(dim(FCMAT2.FCoverSE))
  file <- paste0(dir,"FCMAT2.FCoverSE.",dataset,".RData")
  save(FCMAT2.FCoverSE,file=file)

  FCMAT2.PV <- matpv
  print(dim(FCMAT2.PV))
  file <- paste0(dir,"FCMAT2.PV.",dataset,".RData")
  save(FCMAT2.PV,file=file)

  FCMAT2.SE <- matse
  print(dim(FCMAT2.SE))
  file <- paste0(dir,"FCMAT2.SE.",dataset,".RData")
  save(FCMAT2.SE,file=file)
}
