#' Map reference chemicals to super targets
#'
#' @param sigset Name of the signature set.
#' @param sigcatlog Nmae of the catalog file
#' @return the trimmed signature table
#' @export
mapRefchemToSupertarget <- function(sigset="screen_large",
                                    sigcatalog="signatureDB_master_catalog 2021-02-10",
                                    chem.annotation.file = "../input/chemicals/HTTR_chemical_annotations_2021-02-12.xlsx") {

  printCurrentFunction()

  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  mat <- read.xlsx(file)
  mat = mat[mat[,sigset]==1,]
  chems = read.xlsx(chem.annotation.file)
  names(chems) = tolower(names(chems))
  x = chems$dtxsid
  y = x[duplicated(x)]
  rownames(chems) = chems$dtxsid

  chem.list = unique(chems[,c("dtxsid","casrn","name","use_class","target_all_information","target_gene_family","target_genes")])
  rownames(chem.list) = chem.list$dtxsid

  chem.list$stlist = NA
  st.list = sort(unique(mat$super_target))

  for(i in 1:nrow(chem.list)) {
    dtxsid = chem.list[i,"dtxsid"]
    temp = chems[dtxsid,]
    x1 = temp[1,"target_all_information"]
    x2 = temp[1,"target_gene_family"]
    x3 = temp[1,"target_genes"]
    tlist = str_split(x1,"\\|")[[1]]
    if(!is.na(x2)) tlist = c(tlist,str_split(x2,"\\|")[[1]])
    if(!is.na(x3)) tlist = c(tlist,str_split(x3,"\\|")[[1]])
    tlist = str_replace_all(tlist,":Positive","")
    tlist = str_replace_all(tlist,":Negative","")
    tlist = str_replace_all(tlist," Agonist","")
    tlist = str_replace_all(tlist," Antagonist","")
    tlist = str_replace_all(tlist," Inhibitor","")

    tlist = sort(unique(tlist))
    x = NULL
    for(target in tlist) {
      y = NULL
      if(is.element(target,c("CA"))){
        if(is.element(target,st.list)) y = target
      }
      else y = grep(target,st.list)
      if(length(y)>0) x = c(x,y)
    }
    cat(dtxsid,"\n",tlist,"\n",st.list[x],"\n\n")
    stall = st.list[x]
    stall=paste(stall,collapse="|")
    chem.list[dtxsid,"stlist"]=stall
  }
  cat(nrow(chem.list),"\n")
  x = chem.list$stlist
  x[x==""] = "-"
  chem.list$stlist = x
  chem.list = chem.list[nchar(chem.list$stlist)>0,]
  cat(nrow(chem.list),"\n")
  file = "../input/chemicals/refchem_super_target_map.xlsx"
  write.xlsx(chem.list,file)
}

