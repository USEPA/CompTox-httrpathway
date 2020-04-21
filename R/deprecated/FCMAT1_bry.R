#############################################################################################################
##                                               FCMAT1_bry                                                ##
#############################################################################################################
#' hopefully a faster version of the very for-loopy FCMAT1 found in the myGSEA pagackge taken from DERICK on
#' feb 20 2020.  Written by Bryant Chambers mar 15 2020
#' 
#'about:
#'this script builds the giant FCMAT1 from a folder of tsv files in a directory of your choosing
#'
#'important. it asssumes that you have all files as .tsv s in a single folder. it doesn't check 
#'that all files are tsvs or exclude non tsv files. I just wrote this as a trial and it seems to 
#'significantly outperform the original. (a 3421 file set took 6 min while the original took 8 hours)

###########################################################################################################
#load libraries:
suppressMessages(library(tidyverse))

###########################################################################################################
#load data directories
#!!!!!!!!!!!!!!!!! write the dir to the folder of your input tsvs here: 
filelist <- list.files("~/ipynb/Projects/SMI/SMI_c05GSEA_myGSEArichardbuild/input/myGSEA_input_LINCS_handpicked_allp1andp2/") #list of all files to iterate through
FILE_DIR <- "~/ipynb/Projects/SMI/SMI_c05GSEA_myGSEArichardbuild/input/myGSEA_input_LINCS_handpicked_allp1andp2/" #same as above - should just enter into the function




###########################################################################################################
#read all files and bind them simultaneously 
#this builds a perliminary version of the FCMAT1
allframes <- suppressMessages(lapply(1:length(filelist),function(x)read_tsv(paste0(FILE_DIR,filelist[x]))))
allframesbound <- do.call(rbind, allframes)



###########################################################################################################
#clean up the mat so it look identical to the original fcmat
names(allframesbound) <- c("probe_id","sample_key","basemean","l2fc","se","stat","pvalue","padj")               ## double check these so that they match yours
allframesbound$basemean <- as.integer(allframesbound$basemean)
FCMAT1 <- as.data.frame(allframesbound)


############################################################################################################
#save outputs
save(FCMAT1, file = "FCMAT1_StressModuleEnrichment.RData")