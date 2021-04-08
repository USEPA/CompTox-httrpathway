#' Plot conc response for genes fro one signature and one chemical
#'
#' @param sigset Name of the signature set.
#' @param dataset Name of the data set.
#' @param method Pathway scoring method in c("fc", "gsva", "mygsea")
#' @param bmr_scale	bmr scaling factor. Default = 1.349
#' @param mc.cores Number of cores to parallelize with.
#' @param do.load If TRUE, load the SIGNATURE_CR file, otherwiseassume that it is in memory
#' to.file to.file = T saves the output to a file; otherwise it's returned.
#' @param pval Desired cutoff p-value.
#' @param nametag Optional descriptor tag to attach to file outputs for
#'   experimental/non-default runs.
#' @param plotrange The x-range of the plot as a vector of 2 elements, this can be changed for special cases, but defaults to 0.001 to 100
#'
#' @import data.table
#' @import parallel
#' @import openxlsx
#'
#' @export
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_block_123
#' u2os_toxcast_pfas_pe1_normal
#' PFAS_HepaRG
#' PFAS_U2OS
#' u2os_pilot_pe1_normal_null_pilot_lowconc
#' u2os_toxcast_pfas_pe1_normal_refchems
#' heparg2d_toxcast_pfas_pe1_normal_refchems
#'----------------------------------------------------------------------------------
geneSignatureCheck <- function(to.file=F,
                               do.load=F,
                               dataset="heparg2d_toxcast_pfas_pe1_normal",
                               pval = .05,
                               nametag = "_conthits",
                               dtxsid="DTXSID0020365",
                               signature_parent="LEE_EARLY_T_LYMPHOCYTE",
                               sigcatalog="signatureDB_master_catalog 2021-03-05",
                               sigset="screen_large",
                               plotrange=c(0.001,100)) {

  printCurrentFunction(paste(dataset))

  if(do.load) {
    file <- paste0("../output/gene_conc_resp_summary/GENE_CR_",dataset,"_", pval, nametag ,".RData")
    print(file)
    load(file=file)

    #fix chemical name so it can be part of a file name
    GENE_CR$proper_name = gsub("\\)","",GENE_CR$name)
    GENE_CR$proper_name = gsub("\\(","",GENE_CR$proper_name)
    GENE_CR$proper_name = gsub(":","",GENE_CR$proper_name)
    GENE_CR$proper_name = gsub("%","Percent",GENE_CR$proper_name)
    GENE_CR <<- GENE_CR
    file <- paste0("../input/fcdata/FCMAT2_",dataset,".RData")
    load(file=file)
    FCMAT2 <<- FCMAT2
    file <- paste0("../input/fcdata/CHEM_DICT_",dataset,".RData")
    load(file=file)
    CHEM_DICT <<- CHEM_DICT
  }

  gene_cr = GENE_CR#[is.element(GENE_CR$dtxsid,dtxsid),]
  cname = gene_cr[1,"name"]
  pname = gene_cr[1,"proper_name"]


  for(sid in unique(gene_cr$sample_id)) {
    temp1 = gene_cr[is.element(gene_cr$sample_id,sid),]
    temp2 = CHEM_DICT[is.element(CHEM_DICT$sample_id,sid),]
    sk.list = temp2$sample_key
    temp3 = FCMAT2[sk.list,]
    for(i in 1:nrow(temp1)) {
      gene = temp1[i,"gene"]
      resp2 = temp3[,gene]
      resp1 = temp1[i,"resp"]
      resp1 = str_split(resp1,"\\|")[[1]]
      browser()
    }

  }

  browser()


  if(to.file) {
    fname <- paste0("../output/gene_signature_check/geneSignatureCheck_",dataset,"_",pname,"_",signature_parent,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,2),mar=c(5,5,4,2))

  file = paste0("../input/signatures/",sigcatalog,".xlsx")
  catalog = read.xlsx(file)
  catalog = catalog[catalog[,sigset]==1,]

  file = "../input/signatures/signatureDB_genelists.RData"
  load(file=file)

  signature.list = catalog[is.element(catalog$parent,signature_parent),"signature"]
  for(signature in signature.list) {
    genes = genelists[signature][[1]]
    genes.in = genes
    genes.all = unique(gene_cr$gene)
    genes.out = genes.all[!is.element(genes.all,genes.in)]
    genes.in = genes.in[is.element(genes.in,colnames(FCMAT2))]
    genes.out = genes.out[is.element(genes.out,colnames(FCMAT2))]
    temp = gene_cr[is.element(gene_cr$gene,genes),]
    sid.list = unique(temp$sample_id)
    for(sid in sid.list) {
      cd = CHEM_DICT[is.element(CHEM_DICT$sample_id,sid),]
      sk.list = cd$sample_key
      mat = FCMAT2[sk.list,]
      conc.list = cd$conc
      name.list = c("conc","resp")
      res = as.data.frame(matrix(nrow=length(conc.list),ncol=length(name.list)))
      names(res) = name.list
      res$conc = conc.list
      for(i in 1:length(sk.list)) {
        sk = sk.list[i]
        vals.in = mat[sk,genes.in]
        vals.out = mat[sk,genes.out]
        res[i,"resp"] = mean(vals.in) - mean(vals.out)
      }
      plot(res$resp~res$conc,log="x",xlim=plotrange,ylim=c(-1,1),main=paste(cname,"\n",signature),cex.main=0.9)
      lines(plotrange,c(0,0))

       subframe = temp[is.element(temp$sample_id,sid),]
       subframe$auc = abs(subframe$top) * (3-log10(subframe$bmd))
       subframe = subframe[order(subframe$auc,decreasing=T),]
       for(i in 1:nrow(subframe)){
         geneConcRespPlot(subframe[i,],plotrange=plotrange)
         gene = subframe[i,"gene"]
         if(!to.file) browser()
       }
     }
  }
  if(to.file) dev.off()
}


