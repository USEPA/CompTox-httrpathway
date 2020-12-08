library(tcpl)
#-----------------------------------------------------------------------------------#
#' Get the TOxCast PODs using input from
#-----------------------------------------------------------------------------------#
toxcastPOD <- function(do.prep=F) {
  printCurrentFunction()

  if(do.prep) {
    #tcplConf(user='rjudson', pass='Catman2@', host='mysql-res1.epa.gov', drvr = 'MySQL',db = 'invitrodb')
    #tcplConf(user='rjudson', pass='Catman2@', host='mysql-res1.epa.gov', drvr = 'MySQL',db = 'prod_internal_invitrodb_v3_2')
    tcplConf(user='rjudson', pass='Catman2@', host='ccte-mysql-res.epa.gov', drvr = 'MySQL',db = 'prod_internal_invitrodb_v3_2')
    setDBConn(password="Catman2@")

    #-----------------------------------------------------------------------------------#
    # Extract data from invitrodb
    #-----------------------------------------------------------------------------------#
    cat("1b: Extract data from invitrodb\n")
    file <- "../input/chemicals/httr_chemical_annotations 2020-09-14.xlsx"
    chems <- read.xlsx(file)
    pod.master.casn <- chems$casrn
    pod.master.casn = c(pod.master.casn,"446-72-0","53123-88-9")
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

    length(unique(filter.invitrodb.hitc1.mc5$casn)) #447 casrn

    filter.invitrodb.hitc1.mc5[, hitcsum := sum(hitc), by=list(casn)]
    filter.invitrodb.hitc1.mc5[, modl_ga_uM := ifelse(!is.na(modl_ga), 10^modl_ga, NA)]
    filter.invitrodb.hitc1.mc5[, min_modl_ga_uM := min(modl_ga_uM, na.rm=TRUE), by=list(chid)]
    file='../toxcast/toxcast_master.RData'


    toxcast.master <- as.data.frame(filter.invitrodb.hitc1.mc5)
    save(toxcast.master, file='../toxcast/toxcast_master.RData')

    toxcast.all <- as.data.frame(invitrodb.mc5)
    save(toxcast.all, file='../toxcast/toxcast_all.RData')
  }

  load(file='../toxcast/toxcast_all.RData')
  mat.all = toxcast.all
  load(file='../toxcast/toxcast_master.RData')
  mat <- toxcast.master
  mat0 = mat
  mat <- mat[!is.na(mat$hitpct),]
  chems <- unique(mat[,c("dsstox_substance_id","casn","chnm")])
  names(chems) <- c("dtxsid","casrn","name")
  chems$nhit <- NA
  chems$nassay <- NA
  chems$pod_uM <- NA
  chems$pod_loguM <- NA
  rownames(chems) <- chems$dtxsid
  for(dtxsid in chems$dtxsid) {
    temp <- mat[is.element(mat$dsstox_substance_id,dtxsid),]
    temp0 <- mat.all[is.element(mat.all$dsstox_substance_id,dtxsid),]
    x <- log10(temp$modl_ga_uM)
    q <- quantile(x,probs=seq(0,1,0.05))
    chems[dtxsid,"pod_loguM"] <- q[2]
    chems[dtxsid,"pod_uM"] <- 10**q[2]
    chems[dtxsid,"nhit"] <- length(x)
    chems[dtxsid,"nassay"] <- length(unique(temp0$aenm))
  }
  file <- "../toxcast/toxcast_pod.xlsx"
  write.xlsx(chems,file)
}
