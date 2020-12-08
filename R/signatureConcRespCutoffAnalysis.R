#--------------------------------------------------------------------------------------
#' Run analyses to help determine the hitcall and top_over_cutoff thresholds
#' @param to.file If TRUE, save plots to a file
#' @param do.load If TRUE, load the inpout file to a global
#' @param dataset The name of the data set to use
#' @param sigset THe set of signatures to use
#' @param method signature scoring method in c("fc", "gsva", "mygsea")
#' @param hccut.list The list of hitcall cutoffs to use
#' @param tccut.list The list of top_ver_cutoff cutoffs to use
#' heparg2d_toxcast_pfas_pe1_normal
#' u2os_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123

#--------------------------------------------------------------------------------------
signatureConcRespCutoffAnalysis <- function(to.file=F,
                                            do.load=F,
                                            dataset="mcf7_ph1_pe1_normal_block_123",
                                            sigset="screen_large",
                                            method="fc",
                                            hccut.list=c(0.5,0.6,0.7,0.8,0.9,0.95),
                                            tccut.list=c(1,1.5,2,2.5),
                                            bmd.list=c(100,10,1)) {
  printCurrentFunction(paste(dataset,sigset,method))
  if(to.file) {
    fname <- paste0("../output/signature_conc_resp_cutoffs/scrTrends ",dataset,"_",sigset,"_",method,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(4,4,2,2))

  if(do.load) {
    file = paste0("../output/signature_conc_resp_summary/SIGNATURE_CR_",sigset,"_",dataset,"_",method,"_0.05_conthits.RData")
    print(file)
    load(file=file)
    mat = SIGNATURE_CR
    temp = unique(mat[,c("sample_id","dtxsid")])
    dtxsid.list = temp$dtxsid
    dtemp = dtxsid.list[duplicated(dtxsid.list)]
    mat = mat[is.element(mat$dtxsid,dtemp),]
    MAT <<- mat
  }
  mat = MAT
  mat[is.na(mat$top_over_cutoff),"top_over_cutoff"] = 0
  mat[is.na(mat$bmd),"bmd"] = 1000
  mat[mat$bmd>1000,"bmd"] = 1000
  mat$auc = -log10(mat$bmd/1000) * mat$top_over_cutoff
  name.list = c("dtxsid","name","tccut","hccut","xname","nsig","r2","rmse")
  row = as.data.frame(matrix(nrow=1,ncol=length(name.list)))
  names(row) = name.list
  set1 = NULL
  for(dtxsid in unique(mat$dtxsid)) {
    temp = mat[is.element(mat$dtxsid,dtxsid),]
    sids = unique(temp$sample_id)
    nsid = length(sids)
    name = temp[1,"name"]

    if(nsid>1) {
      for(i1 in 1:(nsid-1)) {
        for(i2 in (i1+1):nsid) {
          sid1 = sids[i1]
          sid2 = sids[i2]
          mat1 = temp[is.element(temp$sample_id,sid1),]
          mat2 = temp[is.element(temp$sample_id,sid2),]
          for(tc in tccut.list) {
            mat1a = mat1[mat1$top_over_cutoff>=tc,]
            for(hc in hccut.list) {
              mat1b = mat1a[mat1a$hitcall>=hc,]
              rownames(mat1b) = mat1b$signature
              rownames(mat2) = mat2$signature
              slist = mat1b$signature
              slist = slist[is.element(slist,mat2$signature)]
              mat1b = mat1b[slist,]
              mat2 = mat2[slist,]
              if(length(slist)>10) {
                x = mat1b[,"auc"]
                y = mat2[,"auc"]
                cat(name,tc,hc,length(slist),"\n")
                res= lm(y~x)
                sr = summary(res)

                row[1,"dtxsid"] = dtxsid
                row[1,"name"] = name
                row[1,"tccut"] = tc
                row[1,"hccut"] = hc
                row[1,"xname"] = paste("TH",tc,hc)
                row[1,"nsig"] = nrow(mat1b)
                row[1,"r2"] = sr$r.squared
                row[1,"rmse"] = sr$sigma
                set1 = rbind(set1,row)
                if(hc==hccut.list[1] && tc==tccut.list[1]) {
                  #main = paste(name,"\n",tc,hc)
                  #plot(y~x,main=main,xlab="AUC(1)",ylab="AUC(2)",xlim=c(0,20),ylim=c(0,20),cex.lab=1.2,cex.axis=1.2,pch=21,cex=0.1)
                  #lines(c(0,100),c(0,100))
                  #text(0,18,paste("R2",format(sr$r.squared,digits=2)),pos=4)
                  #if(!to.file) browser()
                }
              }
            }
          }
        }
      }
    }
  }
  xvals = set1$xname
  yvals = set1$r2
  par(mfrow=c(1,1),mar=c(4,10,2,2))
  boxplot(yvals~xvals,horizontal=T,main=paste(dataset,"\nR2"),las=1,xlab="R2",ylab="",ylim=c(0,1))
  for(x in c(0,0.2,0.4,0.6,0.8,1)) lines(c(x,x),c(0,1000))
  if(!to.file) browser()
  yvals = set1$rmse
  boxplot(yvals~xvals,horizontal=T,main=paste(dataset,"\nRMSE"),las=1,xlab="RMSE",ylab="",ylim=c(0,5))
  for(x in c(0,1,2,3,4,5)) lines(c(x,x),c(0,1000))
  if(!to.file) browser()
  yvals = set1$nsig
  boxplot(yvals~xvals,horizontal=T,main=paste(dataset,"\nActive Signatures"),las=1,xlab="Active Signatures",ylab="",ylim=c(0,3000))
  for(x in c(0,500,1000,1500,2000,2500,3000)) lines(c(x,x),c(0,1000))
  if(!to.file) browser()

  ##################################################################
  mat = MAT
  mat = mat[mat$hitcall>=0.5,]
  mat = mat[mat$top_over_cutoff>=1,]
  sig.list = sort(unique(mat$signature))
  res1 <- unique(mat[,c("sample_id","dtxsid","name")])
  res1 <- res1[order(res1$dtxsid),]
  rownames(res1) <- res1$sample_id
  nchem <- nrow(res1)
  res2 <- matrix(nrow=nchem,ncol=length(sig.list))
  colnames(res2) <- sig.list
  rownames(res2) <- rownames(res1)
  res2[] <- 0

  yvals = NULL
  xvals = NULL
  for(tc in tccut.list) {
    cat("TC:",tc,"\n")
    mat1 = mat[mat$top_over_cutoff>=tc,]
    for(hc in hccut.list) {
      cat("  HC:",hc,"\n")
      mat2 = mat1[mat1$hitcall>=hc,]
       mat3 = mat2
      res2[] <- 0
      reshc = res2
      for(sid in unique(mat3$sample_id)) {
        #cat("    sid:",sid,"\n")
        mat4 = mat3[is.element(mat3$sample_id,sid),]
        slist = unique(mat4$signature)
        res2[sid,slist] = 1
      }
      xval = paste("TH:",tc,hc)
      for(dtxsid in unique(res1$dtxsid)) {
        sids = unique(res1[is.element(res1$dtxsid,dtxsid),"sample_id"])
        nsid = length(sids)
        if(nsid>1) {
          #cat(dtxsid,nsid,"\n")
          #print(sids)
          #if(nsid==1) browser()
          for(i1 in 1:(nsid-1)) {
            for(i2 in (i1+1):nsid) {
              sid1 = sids[i1]
              sid2 = sids[i2]
              r1 = res2[sid1,]
              r2 = res2[sid2,]
              r12 = r1+r2
              n2 = length(r12[r12==2])
              n1 = length(r12[r12>=1])
              if(n1>0)  ratio = n2/n1
              else ratio = 0
              yvals = c(yvals,ratio)
              xvals = c(xvals,xval)
            }
          }
        }
      }
    }
  }
  #par(mfrow=c(1,1),mar=c(4,4,2,2))
  boxplot(yvals~xvals,horizontal=T,main=paste(dataset,"\nHitcall Tanimoto"),las=1,xlab="Tanimoto",ylab="",ylim=c(0,1))
  for(x in c(0,0.2,0.4,0.6,0.8,1)) lines(c(x,x),c(0,1000))

  if(to.file==F) browser()
  else dev.off()
}
