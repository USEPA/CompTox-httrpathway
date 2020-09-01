library(openxlsx)
library(tidyr)
library(stringr)
library(stringi)
library(reshape2)
library(gplots)
library(RColorBrewer)
#--------------------------------------------------------------------------------------
#' Calcualte the stats for the reference chemicals
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' u2os_toxcast_pfas_pe1_normal
#--------------------------------------------------------------------------------------
refchemdbStats <- function(to.file=F,
                           do.load=F,
                           do.newrefchems=F,
                           do.count0=F,
                           do.count1=F,
                           do.plot=F,
                           do.stats=F,
                           dataset="u2os_toxcast_pfas_pe1_normal",
                           sigset="screen_large",
                           method="fc",
                           celltype="U2OS",
                           hccut=0.95) {
  printCurrentFunction(paste(dataset,sigset,method))
  file = "../input/signatures/refchemdb_chem_filtered_unique_with_new.xlsx"
  rcdb = read.xlsx(file)

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    MAT <<- mat
  }
  if(do.newrefchems) {
    mat = MAT
    mat = mat[mat$hitcall>hccut,]
    mat = mat[mat$bmd<1,]
    mat = mat[mat$top_over_cutoff>2,]
    name.list = c("dtxsid","casrn","name","super_target")
    row = as.data.frame(matrix(nrow=1,ncol=length(name.list)))
    names(row) = name.list
    res = NULL
    nmin = 10
    chem.list <- sort(unique(mat$name))
    for(chem in chem.list) {
      temp = mat[is.element(mat$name,chem),]
      st.list = unique(temp$super_target)
      for(st in st.list) {
        temp1 = temp[is.element(temp$super_target,st),]
        if(nrow(temp1)>=nmin) {
          row[1,"dtxsid"] = temp1[1,"dtxsid"]
          row[1,"casrn"] = temp1[1,"casrn"]
          row[1,"name"] = temp1[1,"name"]
          row[1,"super_target"] = temp1[1,"super_target"]
          print(row)
          res = rbind(res,row)
        }
      }
    }
    file = paste0("../output/signature_refchemdb/",celltype,"/new_refchems_",celltype,".xlsx")
    write.xlsx(res,file)
  }
  if(do.count0) {
    cat("get the stats for the whole data set\n")
    mat = MAT
    mat = mat[is.element(mat$super_target,rcdb$super_target),]
    slist = sort(unique(mat$signature))
    nsig = length(slist)
    name.list = c("signature","super_target","n","nhit","nhit10","nhit1")
    res = as.data.frame(matrix(nrow=nsig,ncol=length(name.list)))
    mathit = mat[mat$hitcall>hccut,]
    names(res) = name.list
    for(i in 1:nsig) {
      sig = slist[i]
      cat(sig,"\n")
      res[i,"signature"] = sig

      if(i==1) {
        temp = mat[is.element(mat$signature,sig),]
        n0 = nrow(temp)
      }
      res[i,"n"] = n0
      temp = mathit[is.element(mathit$signature,sig),]
      if(nrow(temp)>1) res[i,"super_target"] = temp[1,"super_target"]
      temp1 = temp[temp$hitcall>0.9,]
      res[i,"nhit"] = nrow(temp1)
      res[i,"nhit10"] = nrow(temp1[temp1$bmd<10,])
      res[i,"nhit1"] = nrow(temp1[temp1$bmd<1,])
      if(i%%100==0) cat("finished:",i," out of ",nsig,"\n")
    }
    file = paste0("../output/signature_refchemdb/",celltype,"/stats_all_",celltype,".xlsx")
    write.xlsx(res,file)
  }

  if(do.count1) {
    cat("get the stats for the rcdb data set\n")
    mat = MAT
    print(nrow(mat))
    mask = paste0(mat$dtxsid,"_",mat$super_target)
    rmask = paste0(rcdb$dtxsid,"_",rcdb$super_target)
    mat = mat[is.element(mask,rmask),]
    print(nrow(mat))

    slist = sort(unique(mat$signature))
    nsig = length(slist)
    name.list = c("signature","super_target","n","nhit","nhit10","nhit1")
    res = as.data.frame(matrix(nrow=nsig,ncol=length(name.list)))
    names(res) = name.list
    for(i in 1:nsig) {
      sig = slist[i]
      cat(sig,"\n")
      temp = mat[is.element(mat$signature,sig),]
      res[i,"signature"] = sig
      res[i,"super_target"] = temp[1,"super_target"]
      res[i,"n"] = nrow(temp)
      temp1 = temp[temp$hitcall>hccut,]
      res[i,"nhit"] = nrow(temp1)
      res[i,"nhit10"] = nrow(temp1[temp1$bmd<10,])
      res[i,"nhit1"] = nrow(temp1[temp1$bmd<1,])
      if(i%%100==0) cat("finished:",i," out of ",nsig,"\n")
    }
    file = paste0("../output/signature_refchemdb/",celltype,"/stats_rcdb_",celltype,".xlsx")
    write.xlsx(res,file)

    MAT1 <<- mat
    file = paste0("../output/signature_refchemdb/",celltype,"/rcdb_signature_cr_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    save(MAT1,file=file)
  }
  if(do.stats) {
    file = paste0("../output/signature_refchemdb/",celltype,"/stats_all_",celltype,".xlsx")
    res0 = read.xlsx(file)
    file = paste0("../output/signature_refchemdb/",celltype,"/stats_rcdb_",celltype,".xlsx")
    res1 = read.xlsx(file)
    rownames(res0) = res0$signature
    rownames(res1) = res1$signature
    slist = unique(res0$signature)
    slist = slist[is.element(slist,res1$signature)]
    res0 = res0[slist,]
    res1 = res1[slist,]
    name.list = c("signature","super_target","hit","hit10","hit1")
    res = as.data.frame(matrix(nrow=nrow(res0),ncol=length(name.list)))
    names(res) = name.list
    for(i in 1:nrow(res0)) {
      sig = res0[i,"signature"]
      st = res0[i,"super_target"]
      res[i,"signature"] = sig
      res[i,"super_target"] = st

      mat = matrix(nrow=2,ncol=2)
      denom0 = res0[i,"n"]
      denom1 = res1[i,"n"]

      num0 = res0[i,"nhit"]
      num1 = res1[i,"nhit"]
      mat[1,1] = num0
      mat[1,2] = denom0 - num0
      mat[2,1] = num1
      mat[2,2] = denom1 - num1
      p = 1
      if(num0>0 && num1>0) {
        x = prop.test(x=mat,alternative="less")
        p = x$p.value
      }
      res[i,"hit"] = p

      num0 = res0[i,"nhit10"]
      num1 = res1[i,"nhit10"]
      mat[1,1] = num0
      mat[1,2] = denom0 - num0
      mat[2,1] = num1
      mat[2,2] = denom1 - num1
      p = 1
      if(num0>0 && num1>0) {
        x = prop.test(x=mat,alternative="less")
        p = x$p.value
      }
      res[i,"hit10"] = p

      num0 = res0[i,"nhit1"]
      num1 = res1[i,"nhit1"]
      mat[1,1] = num0
      mat[1,2] = denom0 - num0
      mat[2,1] = num1
      mat[2,2] = denom1 - num1
      p = 1
      if(num0>0 && num1>0) {
        x = prop.test(x=mat,alternative="less")
        p = x$p.value
      }
      res[i,"hit1"] = p

    }
    file = paste0("../output/signature_refchemdb/",celltype,"/stats_proptest_",celltype,".xlsx")
    write.xlsx(res,file)
  }
  if(do.plot) {
    if(to.file) {
      fname <- paste0("../output/signature_refchemdb/",celltype,"/signature_refchemdb_",celltype,
                      "_",dataset,"_",sigset,"_",method,".pdf")
      pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
    }
    par(mfrow=c(3,2),mar=c(4,4,2,2))
    mat = MAT1
    mat$proper_name = mat$name

    name.list = sort(unique(mat$name))
    for(name in name.list) {
      cat(name,"\n")
      temp1 = mat[is.element(mat$name,name),]
      rst.list = rcdb[is.element(rcdb$name,name),"super_target"]
      st.list <- sort(unique(temp1$super_target))
      st.list <- st.list[is.element(st.list,rst.list)]
      for(st in st.list) {
        cat("  ",st,"\n")
        temp2 <- temp1[is.element(temp1$super_target,st),]
        temp3 <- temp2[temp2$hitcall>hccut,]
        if(nrow(temp3)>0) {
          temp3 = temp3[order(temp3$top,decreasing=T),]
          for(i in 1:nrow(temp3)){
            #cat(i," out of ",nrow(temp3),"\n")
            signatureConcRespPlot(temp3[i,])
            if(!to.file) browser()
          }
        }
      }
    }
    if(!to.file) browser()
    else dev.off()
  }
}
