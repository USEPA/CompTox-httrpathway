#--------------------------------------------------------------------------------------
#' Look for trends in the raw data by plate within plate group
#'
#' @param dir The directory from which to read all of the raw filesatalog file
#--------------------------------------------------------------------------------------
signatureRankScreen.boxplot <- function(to.file=F,
                                        sys.date="2020-05-01",
                                        dataset="DMEM_6hr_screen_normal_pe_1",
                                        sigset="screen_large",
                                        method = "mygsea") {
  printCurrentFunction()

  if(to.file) {
    fname <- paste0("../output/signature_rank/signature_rank_",dataset,"_",sigset,"_",method," boxplot.pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(1,1),mar=c(5,8,4,3))

  file <- paste0("../output/signature_rank/signature_rank_",dataset,"_",sigset,"_",method,".xlsx")
  mat <- read.xlsx(file)

  mat <- mat[!is.na(mat$deseq.first.ontarget.signature.rank),]
  x <- mat$deseq.first.ontarget.signature.class
  y <- mat$deseq.first.ontarget.signature.rank

  x1 <- unique(x)
  y1 <- y[1:length(x1)]
  y1[] <- 0
  for(i in 1:length(x1)) {
    name <- x1[i]
    vals <- y[is.element(x,name)]
    y1[i] <- median(vals)
  }

  x1 <- x1[order(y1)]
  counter <- 0
  for(val in x1) {
    counter <- counter+1
    sc <- as.character(counter)
    if(nchar(counter)==1) sc <- paste0("0",sc)
    x[is.element(x,val)] <- paste0(sc,": ",val)
  }


  boxplot(y~x,horizontal=T,cex.lab=1.2,cex.axis=0.5,
          ylab="",xlab="rank",
          las=2,outcex=0.5,outpch=21,bg="gray",fg="gray",col="gray",ylim=c(0,0.05))
  for(z in c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1))
    lines(c(z,z),c(0,100),col="gray")
  if(to.file) dev.off()
  else browser()
}
