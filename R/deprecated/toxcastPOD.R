#library(cowplot)
#library(data.table)
#library(ggplot2)
#library(ggpubr)
#library(ggthemes)
#library(grid)
#library(gridExtra)
#library(httk) # version 1.8 used
#library(plyr)
#library(RColorBrewer)
#library(reshape2)
#library(RMySQL)
#library(scales)
#library(stringr)
library(tcpl) # version 2.0.1 used
#-----------------------------------------------------------------------------------#
# Get the TOxCast PODs using Katie's code
#
# Katie Paul Friedman, paul-friedman.katie@epa.gov
# Original 27 Mar 2018
# Updated 6 Dec 2018
#-----------------------------------------------------------------------------------#
toxcastPOD <- function(do.step1=F, # needed
                       do.step2=F, # not needed
                       do.step3=F, # not needed
                       do.step4=F, # not needed
                       do.step5=F, # needed
                       do.step6=F, # ?
                       do.step7=F, # ?
                       do.step8=F  # needed
  ) {
  printCurrentFunction()

  # Invitrodb version 3.0, 7Dec18
  tcplConf(user='rjudson', pass='Catman2@', host='mysql-res1.epa.gov', drvr = 'MySQL',db = 'invitrodb')
  tcplConf(user='rjudson', pass='Catman2@', host='mysql-res1.epa.gov', drvr = 'MySQL',db = 'prod_internal_invitrodb_v3_2')
  setDBConn(password="Catman2@")
  #con <- dbConnect(drv = RMySQL::MySQL(), user="",
  #                 password = "",
  #                 host = "mysql-res1.epa.gov", database = invitrodb_v3)

  #getwd()
  #setwd("/APCRA_final_retrospective/R") # set working directory and all others will be relative
  #load(file='../in_vivo/pod_master.RData')

  if(do.step1) {
    #-----------------------------------------------------------------------------------#
    # Load new httk information
    #-----------------------------------------------------------------------------------#
    cat("1a: Load new httk information\n")
    #load('../httk/new-httk-2018-03-26.RData')

    #chem.physical_and_invitro.data <- add_chemtable(new.data2,
    #                                                current.table=chem.physical_and_invitro.data,
    #                                                data.list=list(Compound="Compound",
    #                                                               CAS="CAS",
    #                                                               DSSTox.GSID="DSSTox_Substance_Id",
    #                                                               Clint="Human.Clint",
    #                                                               Clint.pValue="Human.Clint.pValue",
    #                                                               Funbound.plasma="Human.Funbound.plasma",
    #                                                               LogP="logP",MW="MW",
    #                                                               pKa_Accept="pKa_Accept",
    #                                                               pKa_Donor="pKa_Donor"),
    #                                                reference="Unpublished Ceetox",
    #                                                species="Human",
    #                                                overwrite=T)

    #colnames(chem.physical_and_invitro.data)

    #human.httk <- subset(chem.physical_and_invitro.data,
    #                     (!is.na(chem.physical_and_invitro.data[,22])) &
    #                       (!is.na(chem.physical_and_invitro.data[,28]))) #22=='Human.Clint' & 28=='Human.Funbound.plasma'

    #length(unique(human.httk$CAS)) #732 CAS
    #ade.get <-as.data.frame(subset(human.httk, CAS %in% get_cheminfo(species='Human'))) #718 CAS; what gets lost?
    #ade.get <- as.data.table(ade.get)
    #pod.master.casn <- unique(pod.master[,CASRN])

    #-----------------------------------------------------------------------------------#
    # Extract data from invitrodb
    #-----------------------------------------------------------------------------------#
    cat("1b: Extract data from invitrodb\n")
    file <- "../input/chemicals/HTTR Pilot chemical annotation.xlsx"
    chems <- read.xlsx(file)
    pod.master.casn <- chems$casrn
    invitrodb.chem <- tcplLoadChem(field='casn', val=pod.master.casn, exact=TRUE)
    invitrodb.mc5 <- tcplPrepOtpt(tcplLoadData(lvl=5, type='mc', fld='spid', val=invitrodb.chem$spid))
    invitrodb.hitc1.mc5 <- invitrodb.mc5[which(invitrodb.mc5$hitc==1),]

    # do we lose any CASRN in requiring at least one positive hitcall?

    total.casn <-unique(invitrodb.mc5$casn) #448
    toxcast.casn <- unique(invitrodb.hitc1.mc5$casn) #448
    negatives.casn <- setdiff(total.casn, toxcast.casn) #empty

    #-----------------------------------------------------------------------------------#
    # Derive the min AC50 approximation from ToxCast for these chemicals
    #-----------------------------------------------------------------------------------#
    cat("1c: Derive the min AC50 approximation from ToxCast for these chemicals\n")

    # uncertainty mc7 information is needed for hitpercent

    toxcast.m4id.list <- unique(invitrodb.hitc1.mc5[,m4id])

    mc7 <- runQuery("SELECT
                    mc7.*
                    FROM
                    prod_internal_invitrodb_v3.mc7;","prod_internal_invitrodb_v3_2") %>%
      data.table()

    mc7 <- mc7[m4id %in% toxcast.m4id.list]

    # get mc6 caution flag information, add mc6 and mc7 to the hitc1 table

    positive.spids <- invitrodb.hitc1.mc5[,spid]
    invitrodb.mc6 <- tcplPrepOtpt(tcplLoadData(lvl=6, fld='spid', val=positive.spids, type='mc'))
    setDT(invitrodb.mc6)
    mc6_mthds <- invitrodb.mc6[ , .( mc6_mthd_id = paste(mc6_mthd_id, collapse=",")), by = m4id]
    mc6_flags <- invitrodb.mc6[ , .( flag = paste(flag, collapse=";")), by = m4id]

    invitrodb.hitc1.mc5$mc6_flags <- mc6_mthds$mc6_mthd_id[match(invitrodb.hitc1.mc5$m4id, mc6_mthds$m4id)]
    invitrodb.hitc1.mc5[, flag.length := ifelse(!is.na(mc6_flags), count.fields(textConnection(mc6_flags), sep =','), NA)]
    invitrodb.hitc1.mc5$hitpct <- mc7$hit_pct[match(invitrodb.hitc1.mc5$m4id, mc7$m4id)]
    invitrodb.hitc1.mc5$modl_ga_delta <- mc7$modl_ga_delta[match(invitrodb.hitc1.mc5$m4id, mc7$m4id)]

    #-----------------------------------------------------------------------------------#
    # Filtering the ToxCast data using flags and hitpct
    #-----------------------------------------------------------------------------------#
    cat("1d: Filtering the ToxCast data using flags and hitpct\n")
    invitrodb.hitc1.mc5[flag.length < 3, use.me := 1]
    invitrodb.hitc1.mc5[is.na(flag.length), use.me := 1]
    invitrodb.hitc1.mc5[hitpct <= 0.500 & flag.length >= 3, use.me := 0]
    invitrodb.hitc1.mc5[fitc %in% c(36,45), use.me := 0]
    invitrodb.hitc1.mc5[,length.use.curves:= sum(use.me, na.rm=TRUE), by=list(casn)]
    length(unique(invitrodb.hitc1.mc5[length.use.curves==0]$chid))
    invitrodb.hitc1.mc5[length.use.curves==0] #phenobarbital casn 50-06-6 has no usable curves

    filter.invitrodb.hitc1.mc5 <- invitrodb.hitc1.mc5[use.me==1] #down to 46735 curves from 54048 hitc1 curves
    #(46735/54048)*100
    #100-86.4 #lost 13.6% of the curves

    length(unique(filter.invitrodb.hitc1.mc5$casn)) #447 casrn

    filter.invitrodb.hitc1.mc5[, hitcsum := sum(hitc), by=list(casn)]
    filter.invitrodb.hitc1.mc5[, modl_ga_uM := ifelse(!is.na(modl_ga), 10^modl_ga, NA)]
    filter.invitrodb.hitc1.mc5[, min_modl_ga_uM := min(modl_ga_uM, na.rm=TRUE), by=list(chid)]

    # add row for chems that are negative, set to 100 uM (default high)

    negatives.chemdb<-tcplLoadChem(field='casn', val='50-06-6', exact=TRUE)
    negatives.spids<-negatives.chemdb$spid

    negatives.mc5<-tcplPrepOtpt(tcplSubsetChid(tcplLoadData(lvl=5, fld='spid', val=negatives.spids, type='mc')))
    negatives.mc5 <- unique(negatives.mc5[,list(chid,casn,chnm,code,hitc)])
    negatives.mc5 <- negatives.mc5[!hitc==1]
    negatives.casn <- negatives.mc5[,casn]

    toxcast.master <- merge(filter.invitrodb.hitc1.mc5,
                            negatives.mc5,
                            by=c('chid','casn','chnm','code'),
                            all.x=TRUE,
                            all.y=TRUE)

    toxcast.master[casn %in% negatives.casn, length.use.curves := 0]
    toxcast.master[casn %in% negatives.casn, hitcsum := 0]
    toxcast.master[casn %in% negatives.casn, modl_ga := 2]
    toxcast.master[casn %in% negatives.casn, modl_ga_uM := 100]
    toxcast.master[casn %in% negatives.casn, min_modl_ga_uM := 100]
    toxcast.master[,c('hitc.y'):= NULL]

    save(toxcast.master, file='../toxcast/toxcast_master.RData')
  }
  if(do.step2) {
    #-----------------------------------------------------------------------------------#
    # Histogram of hitcall sum
    #-----------------------------------------------------------------------------------#
    cat("2a: Histogram of hitcall sum\n")
    #need to use the hitcsum to understand how many AC50s will be in the distribution

    load(file='../toxcast/toxcast_master.RData')
    hmc5 <- unique(toxcast.master[,list(casn,chnm,hitcsum)])
    count(hmc5, 'hitcsum')
    median(hmc5$hitcsum, na.rm=TRUE) #56
    range(hmc5$hitcsum) #0 1351

    hist.hitcsum <- ggplot(data=hmc5, aes(hmc5$hitcsum))+
      geom_histogram(binwidth=10)+
      theme(panel.border = element_blank(),
            #panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.title = element_text(face='bold', size=20),
            axis.text = element_text(size=18))+
      xlab('Hitcall Sum')+
      ylab('Frequency')+
      scale_x_continuous(breaks=seq(0,1300,200))
    hist.hitcsum

    hmc5[hitcsum < 100]
    hmc5[casn=='50-06-6']
    file.dir <- paste("Figures/", sep="")
    file.name <- paste("/Figure_hitcsum_", Sys.Date(), ".png", sep="")
    file.path <- paste(file.dir, file.name, sep="")
    dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
    png(file.path, width = 4000, height = 4000, res=480)
    hist.hitcsum
    dev.off()
    browser()
  }
  if(do.step3) {
    #-----------------------------------------------------------------------------------#
    # Optional section: diving a little deeper into how to filter invitrodb data
    #-----------------------------------------------------------------------------------#
    cat("3a: Optional section: diving a little deeper into how to filter invitrodb data\n")
    # relationship between hitpct and modl_ga_delta is not clear; correlation but mostly looks like a blob
    # cannot necessarily infer that large modl_ga_delta means a low hitpct

    summary(lm(hitpct ~ modl_ga_delta, data=invitrodb.hitc1.mc5))
    cor.test(invitrodb.hitc1.mc5$hitpct, invitrodb.hitc1.mc5$modl_ga_delta, method='spearman')

    hitpct <- ggplot(invitrodb.hitc1.mc5, aes(x=modl_ga_delta, y=hitpct)) +
      geom_point() +
      theme_few()+geom_smooth(method=lm)
    hitpct

    # how many curves have low hitpct?

    nrow(invitrodb.hitc1.mc5[hitpct>0.625])
    length(unique(invitrodb.hitc1.mc5[hitpct>0.625]$casn))
    nrow(invitrodb.hitc1.mc5[hitpct>0.5])
    count(invitrodb.hitc1.mc5[hitpct<0.5]$mc6_flags)
    # most prevalent patterns are: 17,16,11,6  (259 times) and 17,16,11,7  (130 times)
    count(invitrodb.hitc1.mc5[hitpct<0.3]$mc6_flags)
    nrow(invitrodb.hitc1.mc5[hitpct<0.5 & flag.length > 2])
    range(invitrodb.hitc1.mc5$flag.length, na.rm=T) # 1 to 6 flags

    flags.hitpct <- ggplot(data=toxcast.master, aes(x=factor(flag.length), y=hitpct))+
      geom_boxplot()+
      xlab('Number of Caution Flags on ToxCast Curve Fit')+
      ylab('Hit percent')

    flags.hitpct.unfiltered <- ggplot(data=invitrodb.hitc1.mc5, aes(x=factor(flag.length), y=hitpct))+
      geom_boxplot()+
      xlab('Number of Caution Flags on ToxCast Curve Fit')+
      ylab('Hit percent')

    double.flags.hitpct <- ggarrange(
      flags.hitpct.unfiltered,
      flags.hitpct,
      widths=c(1,1),
      ncol=2, nrow=1,
      align='h',
      labels=c("A","B"),
      font.label=list(size=18, face='bold'))

    file.dir <- paste("Figures/", sep="")
    file.name <- paste("/Figure_double_flags_hitpct_", Sys.Date(), ".png", sep="")
    file.path <- paste(file.dir, file.name, sep="")
    dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
    png(file.path, width = 6000, height = 4000, res=480)
    double.flags.hitpct
    dev.off()
  }
  if(do.step4) {
    #-----------------------------------------------------------------------------------#
    # Merge in A*STAR HIPPTox Data
    #-----------------------------------------------------------------------------------#
    load(file='../other_data/ASTAR/EPA_export_v2.Rdata')

    ### RData file explanation
    # "HIPPTox.mean.export" and "HIPPTox.err.export" are the means and standard errors of the EC10 values.
    # The unit is log10(micromolar).
    # The mean and standard error values were taken from at least 2 biological replicates.

    # "HIPPTox.mean.export" has four columns.
    # The first three columns are the mean EC10 values for the three cell lines, respectively.
    # Please note that NA means either we have less than two biological replicates that passed QC, or the replicates have too big variations (> 10 times difference, i.e. one order of magnitude).
    # If a chemical was found to be "inactive" (i.e., the best fit DRC model is a constantly low model), the EC10 values will be assigned to log10(99999) =~ 4.999996.

    # Finally, the fourth column is the minimum EC10 of the three cell lines, deriving from the following rules:
    # 1) If either two or more cell lines with non-NA values, minimum of the non-NA values;
    # 2) If only one cell line with a "active" non-NA value (i.e. < log10(99999)), the non-NA value;
    # 3) If only one cell line with a "non-active" non-NA values (i.e. = log10(99999)), NA;
    # 4) If all cell lines with NA values, NA.

    # Essentially, NA means the chemical should not be used.

    #-----------------------------------------------------------------------------------#
    # Prepare ASTAR data for comparison with ToxCast data
    #-----------------------------------------------------------------------------------#

    HIPPTox.mean.export$casn <- row.names(HIPPTox.mean.export)
    astar <- as.data.table(HIPPTox.mean.export)
    colnames(astar)
    head(astar)
    astar[,min.ec10.log10um := Minimum]

    astar.rmNA <- astar[!is.na(min.ec10.log10um)]
    astar.rmNA[,min.ec10.um := ifelse(!is.na(min.ec10.log10um), 10^min.ec10.log10um, NA) ]
    nrow(astar.rmNA[min.ec10.um > 99990]) #5 chemicals negative
    astar.rmNA[min.ec10.um >99990, min.ec10.um := NA] # change these to NA now to signify negatives
    astar.set <- as.data.table(astar.rmNA)

    colnames(astar.set)

    astar.set[,c('BEAS2B',
                 'HepG2',
                 'HK2',
                 'Minimum'):=NULL] #remove ancillary columns

    astar.set[,aenm := 'HIPPTox_mean']

    # add more substance identifiers from invitrodb
    chem.db <- tcplLoadChem(field='casn', val=astar.set$casn)
    chem.db <- unique(chem.db[,c('casn', 'chnm', 'chid')]) #62 chemicals
    astar.set$chid <- chem.db$chid[match(astar.set$casn,chem.db$casn)]
    astar.set$chnm <- chem.db$chnm[match(astar.set$casn,chem.db$casn)]

    astar.set <- astar.set[!is.na(min.ec10.um)] # no data to report; reduced set from 62 to 57 substances
    save(astar.set, file='../toxcast/ASTAR_data_6dec2018.RData')
    load('../toxcast/ASTAR_data_6dec2018.RData')
    colnames(astar.set)
  }
  if(do.step5) {
    #-----------------------------------------------------------------------------------#
    # Obtain distribution of AC50 values for each CASRN
    #-----------------------------------------------------------------------------------#
    cat("5: Obtain distribution of AC50 values for each CASRN\n")
    load(file='../toxcast/toxcast_master.RData')
    length(unique(toxcast.master$casn)) #448 confirmed
    ac50.dist <- toxcast.master[, as.list(quantile(modl_ga_uM, probs=c(0,0.05,0.1,0.25,0.5,0.75,0.9,0.95,1),
                                                   na.rm=TRUE)),keyby=casn]

    ac50.dist$hitcsum<-toxcast.master$hitcsum[match(ac50.dist$casn, toxcast.master$casn)]
    ac50.dist$min_modl_ga_uM <- toxcast.master$min_modl_ga_uM[match(ac50.dist$casn, toxcast.master$casn)]
    ac50.dist$chnm <- toxcast.master$chnm[match(ac50.dist$casn, toxcast.master$casn)]

    pctnames <- grep(names(ac50.dist), pattern='%',value=TRUE)
    setnames(ac50.dist,
             pctnames,
             gsub(x=pctnames,
                  pattern='(\\d{1,3})\\%',
                  replacement='AC50p\\1',
                  perl=TRUE))

    ac50.dist[,avg.dist := AC50p5 - min_modl_ga_uM]

    setnames(ac50.dist, 'casn', 'chemcas')
    colnames(ac50.dist)
    setcolorder(ac50.dist, c(1,13,11:12,14,2:10))

    # now merge with astar.set so that we can create AEDs for AC50p5 from ToxCast and ASTAR EC10

    #colnames(astar.set) # 50 substances with data (not NA)

    #invitro.dist <- merge(ac50.dist,
    #                      astar.set[,c('casn','min.ec10.um')],
    #                      by.x=c('chemcas'),
    #                      by.y=c('casn'),
    #                      all.x=TRUE)

    #head(invitro.dist)
    #nrow(invitro.dist[min.ec10.um < AC50p5]) #1 chemical
    #invitro.dist[min.ec10.um < AC50p5] #50-29-3, dichlorophenyltrichloroethane, by about 1 uM
  }
  if(do.step6) {
    #-----------------------------------------------------------------------------------#
    # Plot distribution of AC50s versus 5th-ile for all chemicals
    #-----------------------------------------------------------------------------------#
    cat("6: Plot distribution of AC50s versus 5th-ile for all chemicals\n")
    # make long format data for ggplot2

    invitro.dist.long <- melt.data.table(invitro.dist,
                                         id.vars = c('chemcas', 'chnm','hitcsum','avg.dist'),
                                         measure.vars = c('min_modl_ga_uM', 'AC50p5', 'min.ec10.um'),
                                         variable.name = c('percentile'))

    invitro.dist.long <- invitro.dist.long[order(-avg.dist, chnm)]

    # plot compares ToxCast min AC50, 5th percentile, and the HIPPTox minEC10 (only 50 chem)
    invitro.dist.plot <- ggplot(data=invitro.dist.long,
                                aes(x=value,
                                    y=reorder(strtrim(chnm,30), -avg.dist),
                                    group=factor(percentile))) +
      geom_point(aes(col=factor(percentile)))+
      scale_colour_manual(values=c("#999999","#0072B2", "#E69F00"),
                          labels=c("min ToxCast", "5%", "min HIPPTox")) +
      scale_shape_manual(values=c(19,19, 19),
                         labels=c("min ToxCast", "5%", "min HIPPTox")) +
      scale_x_log10()+
      theme_bw() +
      coord_flip()+
      theme(panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text = element_text(size=10),
            axis.title = element_text(size=12, face='bold'))+
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
      theme(legend.position="bottom", legend.title=element_blank())+
      xlab('log10 micromolar')+
      ylab('Chemical')

    file.dir <- paste("../toxcast/", sep="")
    file.name <- paste("/Figure_min_v_5th", Sys.Date(), ".png", sep="")
    file.path <- paste(file.dir, file.name, sep="")
    dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
    png(file.path, width = 5000, height = 4000, res=480)
    invitro.dist.plot
    dev.off()
  }
  if(do.step7) {
    #---------------------------------------------------------------------#
    # Original hit-rates for chemicals in the case study
    #---------------------------------------------------------------------#
    cat("7: Original hit-rates for chemicals in the case study\n")
    load(file='../toxcast/toxcast_master.RData')
    # calc hit rates

    hit.rates <- invitrodb.mc5[ , list(
      total.assay.screened  = .N, #total number of aeids tested in mc
      active.assay.count  = as.double(lw(hitc==1)),  # active count
      inactive.assay.count  = as.double(lw(hitc==0)),  #inactive count
      active.percent = round((lw(hitc==1)/.N)*100,2), #active percent
      inactive.percent = round((lw(hitc==0)/.N)*100,2) #inactive percent
    ), by = list(chid, chnm, casn)]

    range(hit.rates$total.assay.screened) #211 to 4557
    median(hit.rates$total.assay.screened) #883.5

    #phenobarbital (50-06-6) and nalidixic acid (389-08-2)
    hit.rates[casn %in% c('50-06-6','389-08-2')]
    # phenobarb was screened in 290 assays, positive in 4 that were dropped
    # nalidixic acid was positive in 9 of 798 assays

    range(hit.rates$active.percent, na.rm=TRUE) # 0.16 54.21
    range(hit.rates$active.assay.count, na.rm=TRUE) #1 to 1445
    median(hit.rates$active.assay.count, na.rm=TRUE)
  }

  if(do.step8) {
    #---------------------------------------------------------------------#
    # frequency of assays driving the minimum AC50 by chemical
    #---------------------------------------------------------------------#
    cat("8: frequency of assays driving the minimum AC50 by chemical\n")
    colnames(toxcast.master)
    setkey(toxcast.master, min_modl_ga_uM, modl_ga_uM)
    toxcast.min <- subset(toxcast.master, min_modl_ga_uM == modl_ga_uM) # 457 rows

    freq.aenm <- toxcast.min[,.N, by=aenm]
    freq.aenm <- freq.aenm[!is.na(aenm)] # there were two chemicals that were negative in ToxCast completely that have NA values
    browser()
    colnames(freq.aenm)
    freq.aenm <- freq.aenm[order(-N)]
    freq.aenm[, asid := (paste(tstrsplit(aenm, '_', 1))), by=aenm]
    freq.aenm[asid=='NIS', asid:= 'NHEERL']

    # figure to show
    # create inset table of top frequency aenms by min POD
    hist(freq.aenm$N)
    freq.table <- freq.aenm[N>5]
    colnames(freq.table)
    freq.table <- freq.table[,asid := NULL]

    cols <- matrix(c("#000000"),nrow(freq.table), ncol(freq.table))
    tt <- ttheme_minimal(
      core=list(fg_params = list(col = cols, cex=1.0),
                bg_params = list(col=NA)),
      rowhead=list(bg_params = list(col=NA), fg_params = list(fontface = "bold",
                                                              cex=1.0)),
      colhead=list(bg_params = list(col=NA), fg_params = list(cex=1.0)))
    t <- tableGrob(freq.table, theme = tt)

    freq.aeid<-ggplot(data=freq.aenm, aes(x = reorder(aenm, -N), y = N, fill=asid)) +
      geom_bar(stat = "identity") +
      scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9, 10, 20))+
      theme_bw() +
      theme(
        #panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12, face='bold'))+
      theme(axis.text.y = element_text(family = "sans", face = "bold", size=12),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
      scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000066", "#FF0033", "#336600"),
                        labels=c('ACEA','APR','ATG','BSK','CEETOX','CLD','NCCT','NHEERL','NVS','OT','Tanguay','TOX21'))+
      #scale_fill_brewer(palette='Set3')+
      #scale_fill_hue(name='Assay Source')+
      ylab('Frequency')+
      xlab('Assay Endpoint Id for Min(AC50)')+
      annotation_custom(grob = t, xmin = 80, xmax = 170, ymin = 10, ymax = 15)

    file.dir <- paste("Figures/", sep="")
    file.name <- paste("/Figure_freq_aeid", Sys.Date(), ".png", sep="")
    file.path <- paste(file.dir, file.name, sep="")
    dir.create(path=file.dir, showWarnings = FALSE, recursive = TRUE)
    png(file.path, width = 4000, height = 5000, res=480)
    freq.aeid
    dev.off()

    #---------------------------------------------------------------------#
    # Write the necessary invitro data to RData
    #---------------------------------------------------------------------#

    save(invitrodb.mc5,
         invitrodb.hitc1.mc5,
         filter.invitrodb.hitc1.mc5,
 #        astar.set,
         toxcast.master,
         invitro.dist, file="../toxcast/source_invitro_data.RData")
  }

}
