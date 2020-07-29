#--------------------------------------------------------------------------------------
#' Run teh WGCNA method on FCMAT2 to determine gene sets to use
#'
#' @param dataset The L2fc matrix data set
#' @param dir The directory where the data file lives
#' @param do.read If TRUE, read in FCMAT2 to a gloabal
#' @return
#'
#' heparg2d_toxcast_pfas_pe1_normal
#' mcf7_ph1_pe1_normal_good_pg
#' DMEM_6hr_pilot_normal_pe_1
#--------------------------------------------------------------------------------------
library(WGCNA)
#  ‘impute’, ‘preprocessCore’, ‘GO.db’
FCMAT2.WGCNA <- function(to.file=F,
                         do.read=F,
                         do.wait=T,
                         dataset="heparg2d_toxcast_pfas_pe1_normal",
                         celltype="HepaRG") {
  printCurrentFunction(dataset)
  if(to.file) {
    fname <- paste0("../output/signature_wgcna/wgcna_",dataset,"_",celltype,".pdf")
    pdf(file=fname,width=8,height=10,pointsize=12,bg="white",paper="letter",pagecentre=T)
  }
  par(mfrow=c(3,1),mar=c(4,4,2,2))

  if(do.read) {
    file <- paste0("../input/fcdata/FCMAT2_",dataset,".RData")
    print(file)
    load(file=file)
    FCMAT2 <<-FCMAT2
  }
  cat("file in, start WGCNA\n")
  mat <- FCMAT2
  datExpr <- mat
  #datExpr <- mat[1:100,1:1000]

  #Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  # Plot the results:
  #sizeGrWindow(9, 5)
  #par(mfrow = c(1,2))
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  if(do.wait) browser()

  softPower = 16
  adjacency = adjacency(datExpr, power = softPower);

  cat("step 1\n")
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1-TOM
  geneTree = hclust(as.dist(dissTOM), method = "average");
  plot(geneTree, xlab="", sub="",
       main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04);

  cat("step 2: Module identification using dynamic tree cut\n")
  minModuleSize = 20
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize)
  modules.original <- dynamicMods
  table(dynamicMods)

  cat("step 3: Plot the dendrogram and colors underneath\n")
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  #sizeGrWindow(8,6)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")

  cat("step 4: Calculate eigengenes\n")
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  #sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")

  cat("step 5\n")
  MEDissThres = 0.25
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs

  plotDendroAndColors(geneTree,
                      cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,addGuide = TRUE,
                      guideHang = 0.05)

  cat("step 6: Construct numerical labels corresponding to the colors\n")
  moduleColors = mergedColors
  colorOrder = c("grey", standardColors(50));
  moduleLabels = match(moduleColors, colorOrder)-1
  modules.merged = moduleLabels
  MEs = mergedMEs
  # Save module colors and labels for use in subsequent parts
  file <- paste0("../output/signature_wgcna/wgcna_",dataset,"_",celltype,".RData")
  save(MEs, moduleLabels, moduleColors, geneTree, file = file)

  cat("step 7\n")
  genes <- colnames(FCMAT2)

  modules.original <- paste0("original_",modules.original)
  modules.merged <- paste0("merged_",modules.merged)

  name.list <- c("signature","parent","source","subsource","type","direction","ngene","description","target_class","super_target","color","include0","wgcna")
  nsig <- length(unique(modules.original)) +   length(unique(modules.merged))
  catalog <- as.data.frame(matrix(nrow=nsig,ncol=length(name.list)))
  names(catalog) <- name.list

  sigdb <- list()

  cat("step 8 create the catalog and signature library\n")
  mlabels <- sort(unique(modules.original))
  for(i in 1:length(mlabels)) {
    label <- mlabels[i]
    gi <- genes[modules.original==label]
    sigdb[[i]] <- gi
    sig <- paste0("wgcna_",dataset,"_",celltype,"_",label)
    names(sigdb)[i] <- sig
    catalog[i,"signature"] <- sig
    catalog[i,"parent"] <- sig
    catalog[i,"source"] <- "wgcna"
    catalog[i,"subsource"] <- "wgcna"
    catalog[i,"type"] <- "unidirectional"
    catalog[i,"direction"] <- "unidirectional"
    catalog[i,"ngene"] <- length(gi)
    catalog[i,"description"] <- sig
    catalog[i,"target_class"] <- sig
    catalog[i,"super_target"] <- sig
    catalog[i,"color"] <- "gray"
    catalog[i,"include0"] <- 1
    catalog[i,"wgcna"] <- 1
  }

  mlabels <- sort(unique(modules.merged))
  for(i in 1:length(mlabels)) {
    ii <- i + length(unique(modules.original))
    label <- mlabels[i]
    gi <- genes[modules.merged==label]
    sigdb[[ii]] <- gi
    sig <- paste0("wgcna_",dataset,"_",celltype,"_",label)
    names(sigdb)[ii] <- sig
    catalog[ii,"signature"] <- sig
    catalog[ii,"parent"] <- sig
    catalog[ii,"source"] <- "wgcna"
    catalog[ii,"subsource"] <- "wgcna"
    catalog[ii,"type"] <- "unidirectional"
    catalog[ii,"direction"] <- "unidirectional"
    catalog[ii,"ngene"] <- length(gi)
    catalog[ii,"description"] <- sig
    catalog[ii,"target_class"] <- sig
    catalog[ii,"super_target"] <- sig
    catalog[ii,"color"] <- "gray"
    catalog[ii,"include0"] <- 1
    catalog[ii,"wgcna"] <- 1
  }

  sig <- paste0("wgcna_",dataset,"_",celltype,"_",softPower,"_",minModuleSize)
  file <- paste0("../input/signatures/signatureDB_",sig,".RData")
  save(sigdb,file=file)
  file <- paste0("../input/signatures/signatureDB_",sig,"_catalog.xlsx")
  write.xlsx(catalog,file)

  if(!to.file) browser()
  else dev.off()
}
