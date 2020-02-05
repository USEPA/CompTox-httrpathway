#--------------------------------------------------------------------------------------
#'
#' Build lane plots by chemical list and signature class, across the datasets
#' @param to.file If TRUE, write plots to a file
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#--------------------------------------------------------------------------------------
signatureChemicalDiffAcrossDatasets <- function(to.file=F,
                                                  chemical.target="ER",
                                                  signature.super_class="estrogen",
                                                  pathset="PathwaySet_20191107",
                                                  method = "fc") {
  printCurrentFunction()
  file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
  chems <- read.xlsx(file)
  chems <- chems[is.element(chems$target_key,chemical.target),]
  dtxsid.list <- chems$dtxsid
  nchem <- length(dtxsid.list)

  file <- "../input/processed_signature_data/signature_catalog 2019-11-07.xlsx"
  catalog <- read.xlsx(file)
  catalog <- catalog[catalog$useme==1,]
  catalog <- catalog[is.element(catalog$super_class,signature.super_class),]
  signature.list <- catalog$signature
  signature.list <- sort(signature.list)

  dataset.list = c(
    "DMEM_6hr_pilot_normal_pe_0",
    "DMEM_12hr_pilot_normal_pe_0",
    "DMEM_24hr_pilot_normal_pe_0",

    "DMEM_6hr_pilot_none_pe_0",
    "DMEM_12hr_pilot_none_pe_0",
    "DMEM_24hr_pilot_none_pe_0",

    "DMEM_6hr_pilot_normal_pe_1",
    "DMEM_12hr_pilot_normal_pe_1",
    "DMEM_24hr_pilot_normal_pe_1",

    "DMEM_6hr_pilot_none_pe_1",
    "DMEM_12hr_pilot_none_pe_1",
    "DMEM_24hr_pilot_none_pe_1",

    "PRF_6hr_pilot_normal_pe_0",
    "PRF_12hr_pilot_normal_pe_0",
    "PRF_24hr_pilot_normal_pe_0",

    "PRF_6hr_pilot_none_pe_0",
    "PRF_12hr_pilot_none_pe_0",
    "PRF_24hr_pilot_none_pe_0",

    "PRF_6hr_pilot_normal_pe_1",
    "PRF_12hr_pilot_normal_pe_1",
    "PRF_24hr_pilot_normal_pe_1",

    "PRF_6hr_pilot_none_pe_1",
    "PRF_12hr_pilot_none_pe_1",
    "PRF_24hr_pilot_none_pe_1"
  )
  if(to.file) {
    fname <- paste0("../output/signature_conc_resp_laneplots/signatureChemicalDiffPlotAcrossDatasets_",chemical.target,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(4,2),mar=c(4,4,2,2))

  dataset.list <- dataset.list[1:12]
  #dataset.list <- dataset.list[1:3]

  mat <- NULL

  #signature.list <- signature.list[1:3]
  for(dataset in dataset.list) {
    file <- paste0("../output/signature_conc_resp_summary/PATHWAY_CR_",pathset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    temp <- PATHWAY_CR
    temp <- temp[is.element(temp$dtxsid,dtxsid.list),]
    temp <- temp[is.element(temp$signature,signature.list),]
    temp$dataset <- dataset
    mat <- rbind(mat,temp)
  }
  x <- mat$bmd
  x[is.na(x)] <- 1000
  mat$bmd <- x

  mat$shrinkage <- NA
  mat$pe <- NA
  mat$tim <- NA

  for(dataset in dataset.list) {
    shrinkage <- "none"
    if(contains(dataset,"normal")) shrinkage <- "normal"
    mat[is.element(mat$dataset,dataset),"shrinkage"] <- shrinkage
    pe <- 0
    if(contains(dataset,"pe_1")) pe <- 1
    mat[is.element(mat$dataset,dataset),"pe"] <- pe
    time <- 6
    if(contains(dataset,"12hr")) time <- 12
    if(contains(dataset,"24hr")) time <- 24
  }

  name.list <- c("signature","condition","bmd.i","bmd.j","delta")
  row <- as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(row) <- name.list

   for(dtxsid in dtxsid.list) {
    temp1 <- mat[is.element(mat$dtxsid,dtxsid),]
    name <- temp1[1,"name"]
    cat(name,"\n")
    res <- NULL
    for(signature in signature.list) {
      #cat(signature,"\n")
      temp2 <- temp1[is.element(temp1$signature,signature),c("time","shrinkage","pe","bmd")]

      condition <- "normal-none"
      for(i in 1:nrow(temp2)) {
        for(j in 1:nrow(temp2)) {
          if(temp2[i,"shrinkage"]=="normal" && temp2[j,"shrinkage"]=="none") {
            if(temp2[i,"time"]==temp2[j,"time"] && temp2[i,"pe"]==temp2[j,"pe"]) {
              row[1,"condition"] <- condition
              row[1,"signature"] <- signature
              row[1,"bmd.i"] <- log10(temp2[i,"bmd"])
              row[1,"bmd.j"] <- log10(temp2[j,"bmd"])
              row[1,"delta"] <- log10(temp2[i,"bmd"])-log10(temp2[j,"bmd"])
              res <- rbind(res,row)
            }
          }
        }
      }

      condition <- "pe1-pe0"
      for(i in 1:nrow(temp2)) {
        for(j in 1:nrow(temp2)) {
          #browser()
          if(temp2[i,"pe"]==1 && temp2[j,"pe"]==0) {
            #browser()
            if(temp2[i,"time"]==temp2[j,"time"] && temp2[i,"shrinkage"]==temp2[j,"shrinkage"]) {
              #browser()
              row[1,"condition"] <- condition
              row[1,"signature"] <- signature
              row[1,"bmd.i"] <- log10(temp2[i,"bmd"])
              row[1,"bmd.j"] <- log10(temp2[j,"bmd"])
              row[1,"delta"] <- log10(temp2[i,"bmd"])-log10(temp2[j,"bmd"])
              res <- rbind(res,row)
            }
          }
        }
      }

    }
    res <- unique(res)
    boxplot(res$delta~res$condition,main=name,cex.lab=1.2,cex.axis=1.2)
    lines(c(-100,100),c(0,0),col="black")
    if(!to.file) browser()
  }

  if(to.file) dev.off()
}

