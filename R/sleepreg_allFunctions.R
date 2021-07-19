
# ---------------------------------
#' @title Down-Sample Accelerometer Files
#' @description Down-samples .csv files to 1Hz, increasing speed and creating required input format for 'GGIR_from_csv' function.
#'
#' Minimum required inputs: 'acceldir', 'col.timestamp', 'col.accel'.
#'
#' @param acceldir Directory containing raw .csv accelerometer files
#' @param alloutdir General output directory, default created if not specified
#' @param dsdir Directory for down-sampled files
#' @param col.timestamp Column of raw .csv files containing timestamp
#' @param col.accel Columns of raw .csv files containing x-y-z accelerometer data e.g., c(1:3)
#'
#' @return
#' @export
#'
#' @examples
#' ds_accel_csv(acceldir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/106", col.timestamp = 1, col.accel = c(2,3,4))
#' ds_accel_csv(acceldir = "~/R/SRI_Study/Geneactiv_Data/reducedData", col.timestamp = 1, col.accel = c(2,3,4))
#' ds_accel_csv(acceldir = "~/R/Biobank/Light_Data/265_evenfreq", col.timestamp = 5, col.accel = c(1,2,3))
ds_accel_csv <- function(acceldir = c(),
                         alloutdir = c(),
                         dsdir = c(),
                         col.timestamp = c(),
                         col.accel = c()
                         ){
  # ----------------------------------------------------
  # Define non-def variables & check for / create dir --------
  if (length(acceldir) == 0){
    stop("Error: Specify directory containing accelerometer .csv files")
  } # If acceldir not specified, return error
  if (length(alloutdir) == 0){
    alloutdir <- paste0(acceldir,"_output") # Specify path for all output
  } # If alloutdir not specified, specify based on acceldir
  if (!dir.exists(alloutdir)) {
    dir.create(alloutdir) # Create directory
  } # If no output dir, create one
  if (length(dsdir) == 0){
    dsdir <- paste0(alloutdir, "/ds_output") # Directory for down-sampled files
  } # If dsdir not specified, specify based on acceldir
  if (!dir.exists(dsdir)) {
    dir.create(dsdir)
  } # If no output dir, create one
  if (length(col.timestamp) == 0){
    stop("Error: Specify columns of accelerometer data timestamp")
  } # If col.timestamp not specified, return error
  if (length(col.accel) == 0){
    stop("Error: Specify columns of accelerometer x,y,z, data")
  } # If col.timestamp not specified, return error

  # ----------------------------------------------------
  # Down-sample across file list --------
  fl <- list.files(acceldir) # List files in acceldir
  fl2 <- gsub("_","-",fl); fl2 <- gsub(" ","-",fl2) # Remove underscores
  for (i in 1:length(fl)){ # Loop over files
    wrdir <- paste(dsdir,fl2[i],sep="/")
    if (!file.exists(wrdir)){ # If down-sampled file doesn't already exist
      apptr <- as.data.frame(data.table::fread(paste(acceldir,fl[i],sep="/"))) # Read ppt data
      ts <- apptr[,col.timestamp] # Timestamps
      secl <- floor(mean(tail(ts)) - mean(head(ts))) # time diff. between start and end of recording (s)
      int <- length(ts)/secl # Down-sampling interval
      rind <- round(seq(from=1,to=length(ts),by=int)) # Down-sampling indices
      apptw <- data.frame(t=rep(NA,length(rind)),x=NA,y=NA,z=NA) # Write to new d.f.
      apptw$t <- seq(from=round(ts[1]),by=1,length.out=length(rind))
      apptw[,c(2,3,4)] <- apptr[rind,col.accel]
      data.table::fwrite(apptw,wrdir) # Write to file
      print(paste0("Down-sampled file extracted: ", fl[i]))
    }
  }
  # ----------------------------------------------------
}

# ---------------------------------
#' @title Get Rolling Window Indices
#' @description Fits rolling windows to time series and outputs start/end indices of each window, relative to time series.
#'
#' Note: keep time units consistent between t, window and step
#'
#' @param t Time series vector
#' @param window Rolling window size
#' @param step Step size for sliding window along time series
#'
#' @return
#' @export
#'
#' @examples
#' rollingWindowInd <- function(t=1:10000, window=3600, step=300)
rollingWindowInd <- function(t=c(),
                             window=c(),
                             step=c()
                             ){
  if (length(t)==0) {
    stop("Error: Specify t")
  }
  if (length(window)==0) {
    stop("Error: Specify window")
  }
  if (length(step)==0) {
    stop("Error: Specify step")
  }

  tLen <- t[length(t)]- t[1] # Length of time series

  if (window >= tLen) {
    stop("Error: Window length must not be longer than time series length")
  }
  if (step > window) {
    stop("Error: Step length must not be longer than window length")
  }

  tAbs <- t - t[1] # Time vector starting at zero
  nWin <- ceiling((tLen-window)/step+1) # Number of windows

  stInd <- rep(NA, nWin); enInd <- rep(NA, nWin);
  for (jj in 1:nWin) {

    tSt <- step*(jj-1)
    tEn <- tSt + window
    stInd[jj] <- which.min(abs(tAbs-tSt))

    if (jj == nWin){
      enInd[jj] <- length(t)
    } else {
      enInd[jj] <- which.min(abs(tAbs-tEn))
    }

  }

  if (stInd[length(stInd)] == enInd[length(enInd)]) {
    stInd <- stInd[-length(stInd)]
    enInd <- enInd[-length(enInd)]
  }

  naBl <- (enInd - stInd) < .8*window | (enInd - stInd) > 1.2*window # For cases where windows are too long / short
  stInd[naBl] <- NA; enInd[naBl] <- NA

  returnLi <- list(stInd, enInd)
  return(returnLi)
}

# ---------------------------------
#' @title Non-wear Detection
#' @description Evaluates epoch-by-epoch non-wear status of accelerometer devices using method described in van Hees et. al. (2011). 15 minute epochs are
#' evaluated based on surrounding 60 minute windows (centered at 15min), where standard deviation < 13mg and range < 50mg in at least two accelerometer axes
#' is required for non-wear classification.
#'
#' Minimum required inputs: dsdir, rmc.col.time, rmc.col.acc
#'
#' @param dsdir Directory for down-sampled files
#' @param rmc.col.time Column of down-sampled .csv files containing timestamp
#' @param rmc.col.acc Columns of down-sampled .csv files containing x-y-z accelerometer data e.g., c(1,2,3)
#' @param nwdir Directory for non-wear data, created under parent output directory if not specified
#' @param sdThres Standard deviation threshold for non-wear classification, applied per window
#' @param rngThres Range threshold for non-wear classification, applied per window
#'
#' @return
#' @export
#'
#' @examples nonwear_detect(dsdir = "~/R/Biobank/Light_Data/265_evenfreq_output/ds_output", rmc.col.time = 1, rmc.col.acc = c(2,3,4))
#' @examples nonwear_detect(dsdir = "~/R/Biobank/Light_Data/127_evenfreq", rmc.col.time = 5, rmc.col.acc = c(1,2,3))
nonwear_detect <- function(dsdir=c(),
                           nwdir=c(),
                           rmc.col.time = 1,
                           rmc.col.acc = c(2,3,4),
                           sdThres = 0.12753,
                           rngThres = 0.4905
                           ){
  if (length(dsdir)==0){
    stop("Error: Specify directory of accelerometer data")
  }
  if (length(nwdir)==0){
    nwdir <- paste0(dirname(dsdir), "/nw_output")
  }
  if (!dir.exists(nwdir)){
    dir.create(nwdir)
  }

  if (length(rmc.col.time) == 0){
    stop("Error: Specify columns of accelerometer data timestamp")
  } # If rmc.col.time not specified, return error
  if (length(rmc.col.acc) == 0){
    stop("Error: Specify columns of accelerometer x,y,z, data")
  } # If rmc.col.time not specified, return error

  dsLi <- list.files(dsdir, pattern = ".csv", full.names = TRUE)
  dsLiNames <- list.files(dsdir, pattern = ".csv")
  for (i in 1:length(dsLi)){ # For each down-sampled file
    tryCatch({ # Catch errors in each loop
      wrdir <- paste(nwdir,dsLiNames[i],sep="/")
      if (!file.exists(wrdir)){
        appt <- read.csv(dsLi[i]) # Read in data
        winOnOff <- rollingWindowInd(appt[,rmc.col.time], 3600, 900)
        winOn <- winOnOff[[1]]; winOff <- winOnOff[[2]]

        winOn <- winOn[!is.na(winOn)]; winOff <- winOff[!is.na(winOff)]

        st15NW <- data.frame(ts =((appt[,rmc.col.time][winOff] - appt[,rmc.col.time][winOn])/2 + appt[,rmc.col.time][winOn] - 450),
                             nonwear = FALSE) # Make d.f. of non-wear blocks

        for (j in 1:nrow(st15NW)){ # For each 60 min block, 15 min step
          countAx = 0
          for (k in 1:length(rmc.col.acc)){
            if (countAx < 2){
              winDat <- appt[,rmc.col.acc[k]][winOn[j]:winOff[j]] # Extract data for this window and axis

              sdev <- sd(winDat)
              rng <- abs(range(winDat)[1] - range(winDat)[2])

              if (sdev < sdThres & rng < rngThres){
                countAx <- countAx + 1
              }
            }
          }
          if (countAx >= 2){ # If >=2 axes satisfy criteria, call this block non-wear
            st15NW$nonwear[j] <- TRUE
          }
        }
        write.csv(st15NW, wrdir, row.names = FALSE)  # Write to parent directory
      }
      print(paste0("Non-wear data extracted: ", dsLiNames[i]))

    }, error = function(e) {
      print(e)
    }, finally = {
      next # Skip to next loop iteration
    }) # Other half of the error-catch function
  }
}

# ----------------------------------------------------
#' @title Raster Plot from Sleep-Wake Summary Data
#' @description Extract raster plots from a binary sleep-wake time series (summarized as a reduced-form data frame)
#'
#' @param SWS Sleep-wake summary data (2-column d.f.). Column 1 contains values 1=sleep, 0=wake, NA=NA, end=recording end.
#' Column 2 contains UNIX timestamps (origin=1970-01-01) identifying the start of each associated column 1 value.
#' Final timestamp specifies recording end.
#' @param rasdir Raster output directory, default created if not specified
#' @param pptName Name of participant or file
#' @param tz Timezone (use OlsonNames() for a list of accepted timezone names)
#'
#' @return
#' @export
#'
#' @examples
#' raster_from_SWS(SWS = SWS, rasdir = "~/R/Biobank/SleepRegPackage/diarydata_output/raster_output", pptName = "ppt1")
raster_from_SWS <- function(SWS = c(),
                            rasdir = c(),
                            pptName = c(),
                            tz = "UTC"
                            ){
  # Check and format function inputs --------------
  # Check input params exist
  if (length(SWS) == 0){
    stop("Error: Specify sleep-wake-summary (SWS) data")
  }
  if (length(rasdir) == 0){
    stop("Error: Specify directory for raster output")
  }
  if (length(pptName) == 0){
    stop("Error: Specify participant or file name")
  }

  # Check data input for correct format
  if (ncol(SWS) != 2 | # If not 2 columns
      length(unique(SWS[,1])) > 4 | # If more than 3 unique values in first column
      !all((SWS[,2][2:(length(SWS[,2]))] - SWS[,2][1:(length(SWS[,2])-1)]) > 0) # If second column is not sequential
  ){
    stop(paste0("Error: Incorrect format sleep-wake summary data, ",pptName))
  }

  # Re-label headers
  names(SWS) <- c("trans","t")
  SWS$trans <- as.numeric(SWS$trans)

  # Function body --------------
  maxdays <- ceiling((range(SWS$t)[2] - range(SWS$t)[1])/60/60/24) # Max days of data this ppt will have
  SWS$tmin <- round(SWS$t/60) # Round t to nearest minute
  onind <- which(SWS$trans == 1 | is.na(SWS$trans)) # Index all sleep onset times

  if (onind[length(onind)] == nrow(SWS)){ # If index is the last, remove it
    onind <- onind[-length(onind)]
  }

  slt <- vector() # Define empty vector
  grp <- vector()
  for (i in 1:length(onind)){ # Extract all times (minute intervals) of sleep
    wrt <- SWS$tmin[onind[i]]:SWS$tmin[(onind[i]+1)] # Times to write to vector
    slt <- c(slt, wrt) # Vector of sleep times updated each loop
    if (SWS$trans[onind[i]] == 1 & !is.na(SWS$trans[onind[i]])){
      grp <- c(grp, rep("Sleep", length(wrt))) # Vector of sleep/NA grp updated each loop
    } else {
      grp <- c(grp, rep("NA", length(wrt)))
    }
  }

  sttso <- SWS$tmin[2] - SWS$tmin[1] # Time difference between recording start and first sleep onset

  rdf <- data.frame(grp = grp, t = slt, tabs = (slt - slt[1] + sttso), day = maxdays, tras = (slt - slt[1] + sttso))
  cuts <- seq(from=1440,by=1440,length.out=(maxdays-1))
  for (i in 1:length(cuts)){
    rdf$day[rdf$tabs >= cuts[i]] <- rdf$day[rdf$tabs >= cuts[i]] - 1
    rdf$tras[rdf$tabs >= cuts[i]] <- rdf$tras[rdf$tabs >= cuts[i]] - 1440
  }
  rdf$trash <- rdf$tras/60 # In hours

  as.POSIXct(SWS$tmin[1]*60,origin="1970-01-01",tz=tz)
  rot <- as.POSIXct(SWS$tmin[1]*60,origin="1970-01-01",tz=tz)
  oh <- as.numeric(substr(rot,12,13)) + as.numeric(substr(rot,15,16))/60

  rdf$tH <- rdf$trash + oh

  sq <- seq(from=oh,to=(oh+24),by=4)
  sq[sq >= 24] <- sq[sq >= 24] - 24
  abb_x <- sq
  abb_y <- rev(seq(from=max(rdf$day),to=1,by=-1))

  rdf$grp <- factor(rdf$grp, levels = c("Sleep", "NA"))

  rst <- ggplot(rdf,aes(x=tH, y=day, color = grp)) +
    geom_point(size=6, shape="\u007C") +
    scale_y_continuous(breaks=seq(from=max(rdf$day),to=1,by=-1),labels=abb_y) +
    scale_x_continuous(breaks=seq(from=oh,to=(oh+24),by=4),limits=c(oh,(oh+24)),labels=abb_x) +
    scale_color_manual(values=c("#000000", "#FFABAB")) +
    xlab(paste0("Time (",tz,", h)")) +
    ylab("Day") +
    theme_classic() +
    theme(text = element_text(size = 18),
          axis.text = element_text(size = 16),
          legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(shape = 15)))

  w <- 200 # Width in mm
  h <- 7/12*w*(1+5/42*(maxdays-7)) # Height (depends on number of days of data)

  ggsave(paste0(rasdir,"/",pptName,".jpg"), device='jpeg',
         plot=rst, width=w, height=h, units="mm", dpi=1400,
         limitsize = FALSE)
  print(paste0("Raster plot saved: ", pptName))
}

# ---------------------------------
#' @title Apply GGIR to Down-Sampled Accelerometer Files
#' @description Specifies parameters and implements GGIR (Migueles et. al., 2019) across all .csv accelerometer files (frequency = 1Hz)
#' in 'acceldir', extracting sleep-wake predictions and sustained inactivity bouts.
#'
#' Minimum required inputs: 'acceldir'
#'
#' @param dsdir Directory for down-sampled files
#' @param alloutdir General output directory, default created if not specified
#' @param outputdir Directory for GGIR output
#' @param rmc.col.acc Columns of accelerometer data in down-sampled files
#' @param rmc.col.time Column of timestamps in down-sampled files
#' @param tz Timezone (use OlsonNames() for a list of accepted timezone names
#'
#' @return
#' @export
#'
#' @examples
#' GGIR_from_csv(dsdir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/106_output/ds_output")
#' GGIR_from_csv(dsdir = "~/R/Biobank/Light_Data/265_evenfreq_output/ds_output")
#'
GGIR_from_csv <- function(dsdir = c(),
                          alloutdir = c(),
                          outputdir = c(),
                          rmc.col.acc = c(2:4),
                          rmc.col.time = 1,
                          tz = "UTC"
                          ){
  # ----------------------------------------------------
  # Define undefined variables & check for / create dirs --------
  if (length(dsdir) == 0){
    stop("Error: Specify directory containing down-sampled .csv files")
  } # If dsdir not specified, return error
  if (!dir.exists(dsdir)) {
    stop("Error: dsdir doesn't exist")
  } # If no output dir, return error

  if(length(alloutdir) == 0){
    sloca <- unlist(gregexpr("/",dsdir)) # Define alloutdir based on dsdir
    sloc <- sloca[length(sloca)]
    alloutdir <- substr(dsdir,1,(sloc-1))
  }

  if (length(outputdir) == 0){
    outputdir <- paste0(alloutdir, "/GGIR_output") # Dir of GGIR output folders
  } # If outputdir not specified, define based on alloutdir
  if (!dir.exists(outputdir)) {
    dir.create(outputdir) # Create directory
  } # If no output dir, create one

  # ----------------------------------------------------
  # Run GGIR across downsampled files --------
  file_list <- list.files(dsdir, pattern = "*.csv", full.names = TRUE) # List of downsampled files
  study_list <- list.files(dsdir, pattern = "*.csv", full.names = FALSE)
  for (k in 1:length(file_list)){
    tryCatch({ # Start of code to catch any error in a loop iteration, write error to SRI file, and skip to next loop iteration
      # --------------------
      # [f] Define 'studyname' -------
      studyname <- study_list[k]
      # --------------------
      # [nf] Check if outputdir already exists -------
      checkdir <- paste(outputdir,"/output_",studyname,sep="")
      if (dir.exists(checkdir)) { # If the output folder directory already exists, skip participant
        next
      }

      # ----------
      # Specify n rows of accel. data
      rmc.nrow <- nrow(data.table::fread(file_list[k])) # Read in data, count number of rows

      # ---------------------------
      # [nf] Run GGIR --------
      ## At some point, go through and delete all the parameters here that are defaults (minor issue)
      ## Are we only specifying that GGIR look at one column of data? -> if so, we can improve the
      ## speed of the down-sampling function by requiring it to only output one accel column and one t col
      GGIR::g.shell.GGIR(
        # General Parameters
        mode=c(1,2,3,4),
        datadir=file_list[k],
        outputdir=outputdir,
        f0=1,
        f1=c(),
        studyname=studyname,
        overwrite = FALSE, # If you want to overwrite the previous milestone data
        do.imp=TRUE, # Do imputation (recommended)
        idloc=2, # id location (1 = file header, 2 = filename)
        storefolderstructure=TRUE,
        chunksize = 1,
        dynrange = c(8),
        # minimumFileSizeMB = 200, # minimum size for analysis
        # desiredtz = "Australia/Melbourne",

        # CSV Read Parameters
        rmc.file = file_list[k], # Filename of file to be read.
        rmc.nrow = rmc.nrow, # Number of rows to read, same as nrow argument in read.csv and in fread.
        rmc.skip = 1, # Number of rows to skip, same as skip argument in read.csv and in fread.
        rmc.dec=".", # Decimal used for numbers, same as skip argument in read.csv and in fread.
        rmc.firstrow.acc = 2, # First row (number) of the acceleration data.
        rmc.firstrow.header= c(), # First row (number) of the header. Leave blank if the file does not have a header.
        rmc.header.length = c(), # If file has header, specify header length (numeric).
        rmc.col.acc = rmc.col.acc, # Vector with three column (numbers) in which the acceleration signals are stored
        rmc.col.temp = c(), # Scalar with column (number) in which the temperature is stored.
        # Leave in default setting if no temperature is avaible. The temperature will be used by g.calibrate.
        rmc.col.time = rmc.col.time, # Scalar with column (number) in which the timestamps are stored. Leave in default setting if timestamps are not stored.
        rmc.unit.acc = "g", # Character with unit of acceleration values: "g", "mg", or "bit"
        rmc.unit.temp = "C", # Character with unit of temperature values: (K)elvin, (C)elsius, or (F)ahrenheit
        rmc.unit.time = "UNIXsec", # Character with unit of timestamps: "POSIX", "UNIXsec" (seconds since origin, see argument origin),
        # "character", or "ActivPAL" (exotic timestamp format only used in the ActivPAL activity monitor).
        rmc.format.time = "%Y-%m-%d %H:%M:%OS", # Format of timestamp, only used for rmc.unit.time: character and POSIX.
        rmc.bitrate = c(), # Numeric: If unit of acceleration is a bit then provide bit rate, e.g. 12 bit.
        rmc.dynamic_range = c(), # Numeric, if unit of acceleration is a bit then provide dynamic range deviation in g from zero,
        # e.g. +/-6g would mean this argument needs to be 6.
        rmc.unsignedbit = FALSE, # Boolean, if unsignedbit = TRUE means that bits are only positive numbers.
        # if unsignedbit = FALSE then bits are both positive and negative.
        rmc.origin = "1970-01-01", # Origin of time when unit of time is UNIXsec, e.g. 1970-1-1
        rmc.desiredtz = tz, # Timezone in which device was configured and experiments took place. If experiments took place in a
        # different timezone, then use this argument for the timezone in whcih the experiments took place and
        # argument configtz to specify where the device was configured (not implemented yet).
        rmc.sf = 1, # Sample rate in Hertz, if this is stored in the file header then the that will be used instead.
        rmc.headername.sf = c(), # If file has a header: Row name (character) under which the sample frequency can be found.
        rmc.headername.sn = c(), # If file has a header: Row name (character) under which the serial number can be found.
        rmc.headername.recordingid = c(), # If file has a header: Row name (character) under which the recording ID can be found.
        rmc.header.structure = c(), # Character used to split the header name from the header value, e.g. ":" or " "
        rmc.check4timegaps = FALSE, # Boolean to indicate whether gaps in time should be imputed with zeros. Some sensing equipment provides
        # accelerometer with gaps in time. The rest of GGIR is not designed for this, by setting this argument to TRUE the
        # gaps in time will be filled with zeros.
        rmc.col.wear = c(), # If external wear detection outcome is stored as part of the data then this can be used by GGIR. This argument specifies
        # the column in which the wear detection (Boolean) is stored.
        rmc.doresample=FALSE, # Boolean to indicate whether to resample the data based on the available timestamps and extracted sample rate from the
        # file header

        #PART 1 Parameters
        # Key functions = reading file, auto-calibration and extracting features
        dayborder= 0,
        windowsizes=c(15,900,3600), #(Epoch length, non-wear detection resolution, non-wear detection evaluation window
        do.cal=TRUE, #Apply autocalibration (recommended)
        do.enmo=TRUE, #Need for physical activity analysis
        do.anglez=TRUE, #needed for sleep detection
        printsummary=FALSE,

        #PART 2 Parameters
        #Key Functions = Non-wear detection, imputation, and basic descriptives
        strategy=1,#If strategy is set to value 1, then check out arguments hrs.del.start and hrs.del.end.
        #If strategy is set to value 3, then check out arguments ndayswindow.
        hrs.del.start=0, #Only relevant when strategy = 1. How many HOURS need to be ignored at the START of the measurement
        hrs.del.end = 0, #Only relevant when strategy = 1. How many HOURS need to be ignored at the END of the measurement
        maxdur=30, #How many Days of measurement do you maximally expect?
        includedaycrit=16, #Number of minimum valid hours in a day to attempt physical activity analysis
        M5L5res=10, #resolution in minutes of M5 and L5 calculation
        winhr=c(5,5), #Size of M5 and L5 (5 hours by default)
        qlevels= c(), #Quantiles to calculate, set value at c() if you dont want quantiles
        qwindow=c(0,24), #Window over which to calculate quantiles
        ilevels=c(), #c(0,100,400,8000), # Acceleration values (metric ENMO) from which a frequency distribution need to be
        iglevels=TRUE,
        mvpathreshold=c(100), # moderate and vigorous physical activity threshold
        bout.metric=4,
        do.parallel = TRUE,

        #PART 3 Parameters
        #Key functions = Sleep detection
        anglethreshold=5,
        timethreshold=5,
        ignorenonwear=TRUE, # If TRUE, non-wear is not detected as sleep (if FALSE, then it will work with imputed data)
        do.part3.pdf=FALSE,# make false inn final product.

        #PART 4 Parameters
        #Key Functions = integrating sleep log (if have one) with sleep detection, storing day and person specific summaries
        excludefirstlast=FALSE, #Exclude first and last night for sleep analysis
        includenightcrit=16, # Number of minimum valid hours in a day to attempt sleep analysis
        def.noc.sleep=c(1), # method as defined by van hees (2018)
        do.visual = FALSE, #not for final analysis
        #  outliers.only = FALSE,

        #Report Generation
        do.report=c(2,4), # What parts does a report need to be generated for (options 2, 4 and 5 )
        visualreport = FALSE,
        dofirstpage = FALSE, # First page of the pdf with simple summary histograms
        viewingwindow = 2 # viewing window of the visual report (1 centres at day and 2 centres at night)
      )
    }, error = function(e) {
      print(e)
    }, finally = {
      next # Skip to next loop iteration
    }) # Other half of the error-catch function
  }
  # ----------------------------------------------------
}

# ---------------------------------
#' @title Calculate Sleep Regularity Index (SRI) from GGIR Output
#' @description Uses sleep windows and sustained inactivity bouts from GGIR output to calculate Sleep Regularity Index scores. Accounts for
#' multiphasic and broken sleep by identifying periods of 'wake' during GGIR-defined sleep windows and periods of 'napping' outside
#' GGIR-defined sleep windows. Uses sustained inactivity bouts to exclude days where sleep onset and offset times are likely
#' miscalculated. Runs across all "output_xxx" directories within 'outputdir', accounting for both multi-file and single-file GGIR output
#' structures.
#'
#' Additional outputs: sleep-wake raster plots, sleep variables, summary of miscalculated nights, sleep-wake vectors.
#'
#' Minimum required inputs: 'outputdir', 'nwdir' (if use.customnonwear = TRUE (default))
#'
#' @param outputdir Directory of GGIR output
#' @param alloutdir General output directory, default created if not specified
#' @param nwdir Directory of non-wear data - specify if using custom (i.e., non-GGIR) non-wear detection. GGIR non-wear data is used preferentially, if available
#' @param use.naps Specify whether 'naps' are included in SRI calculation
#' @param use.WASO Specify whether 'wake after sleep onset' (WASO) periods are included in SRI calculation
#' @param use.miscal Specify whether to filter out nights of data likely miscalculated by GGIR
#' @param use.GGIRnonwear Specify whether to use GGIR's inbuilt non-wear detection
#' @param use.customnonwear Specify whether to use custom (i.e., non-GGIR) non-wear detection (based on van Hees et. al., 2011)
#' @param nonWearInGGIRsleep Specify whether non-wear periods within GGIR's 'sleep windows' are to be included as non-wear
#' @param wr.SWV Specify whether Sleep-Wake Vectors (SWV) are output to file
#' @param wr.raster Specify whether sleep-wake raster plots are output to file
#' @param minSRIdays Minimum number of days of overlapping data to calculate valid SRI scores
#'
#' @return
#' @export
#'
#' @examples
#' SRI_from_GGIR(outputdir = "~/R/Biobank/Light_Data/127_output/GGIR_output")
#' SRI_from_GGIR(outputdir = "~/R/Biobank/GGIR_Version_Testing/105_evensample_output/GGIR_output",
#' nwdir = "~/R/Biobank/GGIR_Version_Testing/105_evensample_output/nw_output")
#' SRI_from_GGIR(outputdir = "~/R/Biobank/Light_Data/127_evenfreq_output/GGIR_output",
#' nwdir = "~/R/Biobank/Light_Data/127_evenfreq_output/nw_output")
SRI_from_GGIR <- function(outputdir = c(),
                          alloutdir = c(),
                          nwdir = c(),
                          use.naps = TRUE,
                          use.WASO = TRUE,
                          use.miscal = TRUE,
                          use.GGIRnonwear = TRUE,
                          use.customnonwear = TRUE,
                          nonWearInGGIRsleep = FALSE,
                          wr.SWV = TRUE,
                          wr.raster = TRUE,
                          minSRIdays = 5
                          ){
  # ----------------------------------------------------
  # Define undefined variables & check for / create dirs --------
  if (length(outputdir) == 0){
    stop("Error: Specify directory containing GGIR output")
  }
  if (!dir.exists(outputdir)) {
    stop("Error: outputdir doesn't exist")
  }

  if(length(alloutdir) == 0){
    sloca <- unlist(gregexpr("/",outputdir)) # Define alloutdir based on outputdir
    sloc <- sloca[length(sloca)]
    alloutdir <- substr(outputdir,1,(sloc-1))
  }

  if (use.customnonwear & length(nwdir) == 0){
    stop("Specify directory containing non-wear data (generated by 'nonwear_detect'). Otherwise, set use.customnonwear = FALSE to rely on GGIR-calculated non-wear only")
  }

  # ----------------------------------------------------
  # Set up SRI.csv file to write to --------
  SRIfile <- paste(alloutdir,"SRI.csv",sep="/")
  SRIheader <- c("file","SRI","SRI_days","SRI_pctl","miscalnights")
  if (!file.exists(SRIfile)) {
    write.table(t(SRIheader),SRIfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
  }

  # ----------------------------------------------------
  # Set up file for miscalculated nights --------
  misclheader <- c("file", "night","preon","poston","preoff","postoff")
  misclfile <- paste(alloutdir,"miscalculated_nights.csv",sep="/")
  if (!file.exists(misclfile)) {
    write.table(t(misclheader),misclfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
  }

  # ----------------------------------------------------
  # Set up file to write Sleep-Wake Vector (SWV) --------
  SWVheader <- c("file","SRItype")
  SWVfile <- paste(alloutdir,"SWV.csv",sep="/")

  if (wr.SWV) {
    if (!file.exists(SWVfile)) {
      write.table(t(SWVheader),SWVfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
    }
  }

  # ----------------------------------------------------
  # Set up dir for raster plots --------
  if (wr.raster == TRUE){
    rasdir <- paste0(alloutdir,"/raster_output")
    if (!dir.exists(rasdir)){
      dir.create(rasdir) # Create directory
    }
  }

  # -----------
  # List GGIR folders -----------
  dir_list <- list.dirs(outputdir, recursive=FALSE)
  dir_list <- dir_list[grep("output_", dir_list)]
  if(length(dir_list) == 0){
    stop("Error: No GGIR output directories available")
  }
  study_list <- list.dirs(outputdir, recursive=FALSE, full.names = FALSE)
  study_list <- study_list[grep("output_", study_list)]

  # -----------
  # Define vector to fill row for error cases ----------
  na_vec <- rep(NA,(length(SRIheader)-1))

  # ----------------------------------------------------
  # Parameters -------
  ## Eventually add these to the list of function input params, but need to place restrictions on range
  WASOmin <- 30
  napmin <- 30
  stepsize <- 2 # step size for rolling window (must be a factor of 1440)
  perc <- 0.95 # percentage of window required to be SIB
  misclwindvec <- c(1.5) # windows to check (in hours), applied for all four pre-/post- onset/offset windows. Outputs SRI scores for last choice only.
  misclperc_outsl <- 0.85
  misclperc_insl <- 0.2
  exclNAhrs <- 6 # Hours of NA (i.e., non-wear) data per day above which the entire day will be excluded

  # load("data/quant.RData")

  # ----------------------------------------------------
  # Loop over files, extracting SWVs and SRI -----------
  for (k in 1:length(dir_list)){
    # Determine unique participants in this iteration of dir_list. For non-parallel GGIR, will only be one. For parallel GGIR, will be multiple
    tryCatch({ # Start of code to catch any error in a loop iteration, write error to SRI file, and skip to next loop iteration

      studyname <- sub("output_","",study_list[k])
      outputfolder <- dir_list[k] # Output folder directory
      nightsumloc <- paste(outputfolder,"/results/QC/part4_nightsummary_sleep_full.csv",sep="") # Nightsummary location
      temp <- read.csv(nightsumloc, header = TRUE)
      nightsummary.ppts <- temp$filename
      temp <- c()
      ppts <- unique(nightsummary.ppts) # List of participants within this GGIR output folder

    for (p in 1:length(ppts)){
      tryCatch({ # Start of code to catch any error in a loop iteration, write error to SRI file, and skip to next loop iteration
        # --------------------
        # Load required data for this ppt -----------
        # Load nightsummary d.f.s
        nightsummary.rd <- read.csv(nightsumloc, header = TRUE)[,c(1,2,3,4,19,20,26,27)] # Read night summary data
        nightsummary <- nightsummary.rd[nightsummary.rd$filename == ppts[p],-1]

        nightsummary2.rd <- read.csv(nightsumloc, header = TRUE)[,c(1,2,3,4,5,7,8,9,14,26,27)] # Different subset of nightsummary for sl. var. extraction
        nightsummary2 <- nightsummary2.rd[nightsummary2.rd$filename == ppts[p],-1]

        # Load SIBs
        sibdir <- paste0(outputfolder,"/meta/ms3.out")
        sibdir.fl <- list.files(sibdir)
        if (length(ppts) == 1){
          sibloc <- paste(sibdir,sibdir.fl[grepl(studyname,sibdir.fl)],sep="/") # For series GGIR output
          load(sibloc) # Load sustained inactivity data to the GE as "sib.cla.sum"
        } else {
          sibloc <- paste(sibdir,sibdir.fl[grepl(ppts[p],sibdir.fl)],sep="/") # For parallel output, find file with ppt name inside
          load(sibloc)
        }

        # Load meta
        metadir <- list.files(outputfolder, recursive = TRUE, pattern = c("meta_", "RData", studyname), full.names = TRUE)
        pptMetadir <- metadir[grepl(ppts[p], metadir)]
        load(pptMetadir)

        # Load config file
        confloc <- list.files(outputfolder, pattern = "config", full.names = TRUE)
        conf <- read.csv(confloc)

        # Timezone
        tz <- conf$value[conf$argument == "desiredtz"]

        # -----------------------------------------------------------
        # Sleep-wake vector calc
        # --------------------
        # -> Get first and last timestamps in minutes ---------
        originIso <- parsedate::parse_iso_8601(M$metashort$timestamp[1])
        originmin <- as.numeric(originIso)/60
        endmin <- as.numeric(parsedate::parse_iso_8601(M$metashort$timestamp[length(M$metashort$timestamp)]))/60
        originclock <- as.numeric(substr(originIso,12,13)) + (as.numeric(substr(originIso,15,16)))/60 + (as.numeric(substr(originIso,18,19)))/60/60 # 24h clock time of first timestamp
        maxdays <- ceiling((endmin-originmin)/60/24)
        wakeonly_vector <- rep(0,ceiling(endmin - originmin)) # Define wake-only vector in 1min increments between recording start and end

        # -> Define SWV based on GGIR sleep times only ---------
        dateSplit <- strsplit(nightsummary$calendar_date, "/") # Extract and re-format date information
        dateSplitDf <- as.data.frame(do.call("rbind", dateSplit))
        dates <- paste(dateSplitDf[,3], dateSplitDf[,2], dateSplitDf[,1], sep="-")

        rftsOn <- rep(NA, length(dates))
        rftsOff <-rep(NA, length(dates))
        for (i in 1:length(dates)) {
          if (nightsummary$sleeponset[i]<24) {
            rftsOn[i] <- paste(dates[i], nightsummary$sleeponset_ts[i], sep=" ")
          } else {
            newdate <- substr((as.POSIXlt.character(dates[i]) + 60*60*26),1,10)
            rftsOn[i] <- paste(newdate, nightsummary$sleeponset_ts[i], sep=" ")
          }
          if (nightsummary$wakeup[i]<24) {
            rftsOff[i] <- paste(dates[i], nightsummary$wakeup_ts[i], sep=" ")
          } else {
            newdate <- substr((as.POSIXlt.character(dates[i]) + 60*60*26),1,10)
            rftsOff[i] <- paste(newdate, nightsummary$wakeup_ts[i], sep=" ")
          }
        }

        isoTsOn <- GGIR::chartime2iso8601(rftsOn, tz=tz) # Re-format times as iso8601
        isoTsOff <- GGIR::chartime2iso8601(rftsOff, tz=tz)

        nightsummary$onsetmin <- round(as.numeric(parsedate::parse_iso_8601(isoTsOn))/60) - originmin
        nightsummary$offsetmin <- round(as.numeric(parsedate::parse_iso_8601(isoTsOff))/60) - originmin

        nightsummary <- nightsummary[(nightsummary$offsetmin - nightsummary$onsetmin) < 16*60,] # Exclude all nightly 'sleep' times > 16h
        if (nrow(nightsummary) < 1){ # If nightsummary is now empty, skip to next line and write an error.
          error_vec <- c(ppts[p],na_vec)
          write.table(t(error_vec),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                      append = TRUE, row.names=FALSE) # Write Error to SRI.csv
          next
        }

        onoffvector_cons <- wakeonly_vector
        for (i in 1:nrow(nightsummary)){ # Make all epochs between onset/offset times equal to one (representing sleep)
          onoffvector_cons[(nightsummary$onsetmin[i]):(nightsummary$offsetmin[i])] <- 1
        }

        # -> Write required SIB data to df, re-format  --------------------
        sibonoff <- sib.cla.sum[,c("sib.onset.time","sib.end.time")] # Subset sib onset and offset times
        sibonoff <- sibonoff[!sibonoff$sib.onset.time == "",] # Check for and remove all empty sibonoff values

        siboncont <- round((as.numeric(parsedate::parse_iso_8601(sibonoff[,1])))/60) - originmin # Convert each set to UNIXmin
        siboffcont <- round((as.numeric(parsedate::parse_iso_8601(sibonoff[,2])))/60) - originmin
        sibonoff_min1 <- data.frame(siboncont, siboffcont) # Get a data frame of sib onset/offset in UNIXmin
        sibonoff_min <- sibonoff_min1[sibonoff_min1[,1] > 0,] # Remove any negative value (for SIB onsets outside SRI days)
        sibonoff_min <- sibonoff_min[sibonoff_min[,2] < length(wakeonly_vector),] # Remove all SIB values outside the range of the swv used throughout.

        # -> Define SWV based on GGIR sleep times and WASO (based on SIBs) --------------
        sib_bool <- logical(nrow(sibonoff_min)) # make a vector of boolean zeros
        for (i in 1:nrow(nightsummary)) { # for each day
          for (j in 1:nrow(sibonoff_min)){ # for each row of sib
            if ( # generate a TRUE value for every sibonset/offset pair with either onset or offset within a sleep-wake period
              (
                (sibonoff_min[j,1] >= nightsummary$onsetmin[i])
                &&
                (sibonoff_min[j,1] <= nightsummary$offsetmin[i])
              )
              ||
              (
                (sibonoff_min[j,2] >= nightsummary$onsetmin[i])
                &&
                (sibonoff_min[j,2] <= nightsummary$offsetmin[i])
              )
            ){
              sib_bool[j] <- TRUE # Set equivalent boolean element to TRUE
            }
          }
        }

        sibonoff_min_insl <- sibonoff_min[sib_bool,] # all sets of sib onset/offset within sleep onset/offset times
        sibonoff_min_outsl <- sibonoff_min[!sib_bool,] # all sets of sib onset/offset outside sleep onset/offset times

        sibOFFON_min_insl <- as.data.frame(cbind( # Rearrange vector, for ease of difference calculating. Glue first onset and last offset to top and tail.
          c(nightsummary$onsetmin[1],sibonoff_min_insl[,2]),
          c(sibonoff_min_insl[,1],nightsummary$offsetmin[length(nightsummary$offsetmin)])
        ))
        diff <- sibOFFON_min_insl[,2] - sibOFFON_min_insl[,1] # Calculate the SIB gaps (offset - previous onset)
        acceptedsibWASO <- sibOFFON_min_insl[(diff >= WASOmin & diff < max(nightsummary$wakeup-nightsummary$sleeponset)*60 ),]
        # Keep all SIB gaps greater than WASOmin and smaller than the longest sleep duration (i.e., wake during the day)

        onoffvector_sibs <- onoffvector_cons
        if (nrow(acceptedsibWASO) > 0){
          for (i in 1:nrow(acceptedsibWASO)){ # Make all epochs between accepted SIB gaps equal to zero (representing wake)
            onoffvector_sibs[(acceptedsibWASO[i,1]):(acceptedsibWASO[i,2])] <- 0
          }
        }

        # -> Define two SWVs, adding both naps and naps+WASO to GGIR sleep-wake times -----------
        onoffvector_allSIBoutsl <- wakeonly_vector
        if (nrow(sibonoff_min_outsl) > 0){
          for (i in 1:nrow(sibonoff_min_outsl)){ # Make all epochs within SIBs outside GGIR sleep/wake equal to one (representing sleep)
            onoffvector_allSIBoutsl[(sibonoff_min_outsl[i,1]):(sibonoff_min_outsl[i,2])] <- 1
          }
        }

        iterlength <- (length(onoffvector_allSIBoutsl)-napmin)/stepsize + 1
        napscore <- replicate(iterlength,0) # Define an empty wake-only vector (all zeros)
        for (i in 1:iterlength){
          if ( sum(onoffvector_allSIBoutsl[((i-1)*stepsize+1):(((i-1)*stepsize+1)+napmin-1)])/napmin >= perc ){
            napscore[i] <- 1
          }
        }

        naploc <- which(napscore %in% 1) # gives vector location of all naps

        onoffvector_GGIRnapwind <- onoffvector_cons
        onoffvector_GGIRWASOnapwind <- onoffvector_sibs
        if (length(naploc) > 0){
          for (i in naploc){
            napepochs <- ((i-1)*stepsize+1):(((i-1)*stepsize+1)+napmin)
            onoffvector_GGIRnapwind[napepochs] <- 1
            onoffvector_GGIRWASOnapwind[napepochs] <- 1
          }
        }

        # -----------------------------------------------------------
        # Sleep variables
        # --------------------
        # -> Extract sleep vars (nightsummary,TST_custom,SE_GGIR,SE_Custom,WASO_GGIR,WASO_Custom,naps)) ---------
        nightsummary2 <- nightsummary2[nightsummary2$night %in% nightsummary$night,] # Exclude all nights not in 'nightsummary'

        SE_GGIR <- nightsummary2[,"SleepDurationInSpt"] / nightsummary2[,"SptDuration"] # Calculate variables directly based on GGIR output
        WASO_GGIR <- nightsummary2[,"SptDuration"] - nightsummary2[,"SleepDurationInSpt"]

        # -> Calculate WASO, TST, SE, based on custom WASO calculation method (>WASOmin) ---------
        # Bugfix ------
        if (nrow(acceptedsibWASO) >= 1){
          sibmeans <- rowMeans(acceptedsibWASO) # which are greater than offset row n and smaller than onset row n + 1
          fix1acceptedsibWASO <- acceptedsibWASO
          for (i in 1:length(sibmeans)){
            if (sum(sibmeans[i] > nightsummary$offsetmin[-nrow(nightsummary)] & sibmeans[i] < nightsummary$onsetmin[-1]) > 0){
              fix1acceptedsibWASO[i,] <- NA
            }
          }
          FIXEDacceptedsibWASO <- fix1acceptedsibWASO[!is.na(fix1acceptedsibWASO[,1]),]
          # ^^ Had to fix a bug with the acceptedsibWASO data.frame - won't affect anything earlier on in the code...
        } else {
          FIXEDacceptedsibWASO <- acceptedsibWASO
        }
        # Extract WASO and naps ------
        WASOdf <- data.frame(night = nightsummary$night,WASO = 0) # create empty data.frame for WASO values
        for (i in 1:nrow(nightsummary)){ # for each night

          WASOoff <- FIXEDacceptedsibWASO[,2][FIXEDacceptedsibWASO[,2] >= nightsummary$onsetmin[i] & FIXEDacceptedsibWASO[,2] <= nightsummary$offsetmin[i]]
          WASOon <- FIXEDacceptedsibWASO[,1][FIXEDacceptedsibWASO[,1] >= nightsummary$onsetmin[i] & FIXEDacceptedsibWASO[,1] <= nightsummary$offsetmin[i]]
          # ^ find all sets of WASO within the on/off times of that night

          if (sum(WASOoff) != 0 & sum(WASOon) != 0){ # If non-zero, write to WASOdf
            WASOdf[i,"WASO"] <- sum(WASOoff-WASOon,na.rm=TRUE)
          }
        }

        WASO_Custom <- (WASOdf[,2])/60
        TST_Custom <- (nightsummary$offsetmin-nightsummary$onsetmin - WASO_Custom)/60
        SE_Custom <- (nightsummary$offsetmin-nightsummary$onsetmin - WASO_Custom) / (nightsummary$offsetmin-nightsummary$onsetmin)

        # Calculate naps
        if (length(naploc) >= 1) {
          endind <- c(which(diff(naploc) != 1),length(naploc)) # Get index of all starts and ends of vector nap locations
          startind <- c(1,which(diff(naploc) != 1)+1)
          naponoff <- data.frame(napstart = rep(NA,length(startind)),napend=NA) # Initialize
          for (i in 1:length(startind)){ # Write nap start and end times to a df.
            naponoff[i,1] <- naploc[startind[i]]
            naponoff[i,2] <- naploc[endind[i]] + napmin
          }
          naps <- rep(0,nrow(nightsummary2)) # Initialize
          for (i in 1:nrow(nightsummary2)){ # Sum all naps within each of the midnight-midnight periods
            naponoffsub <- naponoff[naponoff[,1] < nightsummary2$night[i]*1440 & naponoff[,1] >= nightsummary2$night[i]*1440 - 1440,]
            naps[i] <- (sum(naponoffsub[,2] - naponoffsub[,1]))/60
          }
        } else {
          naps <- rep(0,nrow(nightsummary2))
        }


        # --------------------
        # -> Check for GGIR onset / offset miscalculation --------
        for (i in misclwindvec){
          misclwinddata <- data.frame("participant" = rep(ppts[p],nrow(nightsummary)), night = nightsummary$night,"preon" = NA, "poston" = NA, "preoff" = NA, "postoff" = NA,
                                      "misclwind" = i) # initiate dataframe
          for (j in 1:nrow(misclwinddata)){

            # pre-on
            sibinpreonwind <- sibonoff_min_outsl[(sibonoff_min_outsl$siboffcont > nightsummary$onsetmin[j]-60*i) & (sibonoff_min_outsl$siboffcont <= nightsummary$onsetmin[j]),]
            if (nrow(sibinpreonwind) > 0){ # check if any SIBs are within the misclassification window
              if (sibinpreonwind[1,1] < nightsummary$onsetmin[j]-60*i){ # if first SIBon is outside window, call it start of window
                sibinpreonwind[1,1] <- nightsummary$onsetmin[j]-60*i
              }
              misclwinddata$preon[j] <- sum(sibinpreonwind[,2] - sibinpreonwind[,1])/(60*i)
            } else {
              misclwinddata$preon[j] <- 0
            }

            # post-off
            sibinpostoffwind <- sibonoff_min_outsl[(sibonoff_min_outsl$siboncont < nightsummary$offsetmin[j]+60*i) & (sibonoff_min_outsl$siboncont >= nightsummary$offsetmin[j]),]
            if (nrow(sibinpostoffwind) > 0){ # check if any SIBs are within the misclassification window
              if (sibinpostoffwind[nrow(sibinpostoffwind),2] > nightsummary$offsetmin[j]+60*i){ # if last SIBoff is outside window, call it end of window
                sibinpostoffwind[nrow(sibinpostoffwind),2] <- nightsummary$offsetmin[j]+60*i
              }
              misclwinddata$postoff[j] <- sum(sibinpostoffwind[,2] - sibinpostoffwind[,1])/(60*i)
            } else {
              misclwinddata$postoff[j] <- 0
            }

            sleepmin <- nightsummary$offsetmin[j] - nightsummary$onsetmin[j] # Check that GGIR sleep window is > misclwindow (else leave as NA)
            if (sleepmin > i*60){

              # post-on
              sibinpostonwind <- sibonoff_min_insl[(sibonoff_min_insl$siboncont < nightsummary$onsetmin[j]+60*i) & (sibonoff_min_insl$siboncont >= nightsummary$onsetmin[j]),]
              if (nrow(sibinpostonwind) > 0){ # check if any SIBs are within the misclassification window
                if (sibinpostonwind[nrow(sibinpostonwind),2] > nightsummary$onsetmin[j]+60*i){ # if last SIBoff is outside window, call it end of window
                  sibinpostonwind[nrow(sibinpostonwind),2] <- nightsummary$onsetmin[j]+60*i
                }
                misclwinddata$poston[j] <- sum(sibinpostonwind[,2] - sibinpostonwind[,1])/(60*i)
              } else {
                misclwinddata$poston[j] <- 0
              }

              # pre-off
              sibinpreoffwind <- sibonoff_min_insl[(sibonoff_min_insl$siboffcont > nightsummary$offsetmin[j]-60*i) & (sibonoff_min_insl$siboffcont <= nightsummary$offsetmin[j]),]
              if (nrow(sibinpreoffwind) > 0){
                if (sibinpreoffwind[1,1] < nightsummary$offsetmin[j]-60*i){ # if last SIBoff is outside window, call it start of window
                  sibinpreoffwind[1,1] <- nightsummary$offsetmin[j]-60*i
                }
                misclwinddata$preoff[j] <- sum(sibinpreoffwind[,2] - sibinpreoffwind[,1])/(60*i)
              } else {
                misclwinddata$preoff[j] <- 0
              }
            }
          }
          write.table(misclwinddata[,-7], misclfile, sep = ",", col.names = !file.exists(misclfile), append = TRUE, row.names=FALSE) # write
        }

        misclwinddata[is.na(misclwinddata)] <- 1 # Re-write all NA cases to 1 (cases where window > sleep duration)
        misclnights <- misclwinddata$night[(misclwinddata$preon > misclperc_outsl)|(misclwinddata$poston < misclperc_insl)|
                                             (misclwinddata$preoff < misclperc_insl)|(misclwinddata$postoff > misclperc_outsl)] # Find all miscalculated nights
        miscalculated <- misclwinddata$night %in% misclnights

        # --------------------
        # -> Call missing nights 'NA' in SWVs ---------
        nightSeq <- 1:maxdays
        missingnights <- nightSeq[!(nightSeq %in% nightsummary$night)]
        fMidIn <- round((24-originclock)*60) # Find index of first midnight
        if (length(missingnights) > 0){
          for (i in missingnights){
            NAon <- fMidIn-719+(i-1)*1440
            NAoff <- fMidIn+720+(i-1)*1440

            if (NAoff > length(wakeonly_vector)) { # If offset of last night is longer than SWV length, re-write to SWV length
              NAoff <- length(wakeonly_vector)
            }

            onoffvector_cons[NAon:NAoff] <- NA
            onoffvector_sibs[NAon:NAoff] <- NA
            onoffvector_GGIRnapwind[NAon:NAoff] <- NA
            onoffvector_GGIRWASOnapwind[NAon:NAoff] <- NA
          }
        }

        # --------------------
        # -> Calculate SWVs after removal of miscalculated onset/offset nights --------
        misclnightscount <- length(misclnights) # Number of nights of miscalculated data

        onoffvector_cons_miscl <- onoffvector_cons
        onoffvector_sibs_miscl <- onoffvector_sibs
        onoffvector_GGIRnapwind_miscl <- onoffvector_GGIRnapwind
        onoffvector_GGIRWASOnapwind_miscl <- onoffvector_GGIRWASOnapwind

        if (length(misclnights) > 0){ # Call all 24h periods with miscalculated sleep onset/offset times 'NA'
          for (i in misclnights){

            NAon <- fMidIn-719+(i-1)*1440
            NAoff <- fMidIn+720+(i-1)*1440

            onoffvector_cons_miscl[NAon:NAoff] <- NA
            onoffvector_sibs_miscl[NAon:NAoff] <- NA
            onoffvector_GGIRnapwind_miscl[NAon:NAoff] <- NA
            onoffvector_GGIRWASOnapwind_miscl[NAon:NAoff] <- NA

          }
        }

        # --------------------
        # -> Call non-wear periods NA ---------
        pptTV <- (seq(from=originmin,by=1,length.out=length(onoffvector_cons)))*60 # generate time vector for SWV
        nonWearBool <- rep(FALSE, length(onoffvector_cons)) # Make a boolean non-wear vector, equal length to SWV

        onUNIX <- as.numeric(parsedate::parse_iso_8601(isoTsOn)) # Get sleep window onset and offset times as UNIX
        offUNIX <- as.numeric(parsedate::parse_iso_8601(isoTsOff))
        pptTVoutslBl <- rep(TRUE, length(pptTV))

        if(!nonWearInGGIRsleep){ # If we don't want to include non-wear during GGIR sleep windows:
          for (j in 1:length(onUNIX)) { # Get boolean of all times outside sleep window
            pptTVoutslBl[pptTV > onUNIX[j] & pptTV <= offUNIX[j]] <- FALSE
          }
        }

        if (use.GGIRnonwear) {
          ggirNonWear <- sum(M$metalong$nonwearscore)
          if (ggirNonWear > 0) { # If GGIR non-wear detection has been applied, use these values

            nwTS <- M$metalong$timestamp[M$metalong$nonwearscore >= 2] # Take each TS to represent the centre of a 15 minute window

            nwTSunix <- as.numeric(parsedate::parse_iso_8601(nwTS))

            nwTSunixOn <- nwTSunix - 7.5*60
            nwTSunixOff <- nwTSunix + 7.5*60

            for (j in 1:length(nwTSunix)){
              nonWearBool[pptTV > nwTSunixOn[j] & pptTV <= nwTSunixOff[j]] <- TRUE
            }

            nonWearBool[!pptTVoutslBl] <- FALSE # Make all non-wear inside sleep window = FALSE

            onoffvector_cons[nonWearBool] <- NA # Re-write non-wear times as NA
            onoffvector_sibs[nonWearBool] <- NA
            onoffvector_GGIRnapwind[nonWearBool] <- NA
            onoffvector_GGIRWASOnapwind[nonWearBool] <- NA
            onoffvector_cons_miscl[nonWearBool] <- NA
            onoffvector_sibs_miscl[nonWearBool] <- NA
            onoffvector_GGIRnapwind_miscl[nonWearBool] <- NA
            onoffvector_GGIRWASOnapwind_miscl[nonWearBool] <- NA

          } else { # Else apply custom non-wear detection
            if (use.customnonwear){
              if (length(nwdir) != 0) {
                nwFls <- list.files(nwdir, pattern=".csv", full.names=TRUE) # List all non-wear files

                ppt_redu <- sub(".RData","",ppts[p])
                ppt_redu <- sub(".csv","",ppt_redu)

                appt_nwdir <- nwFls[grepl(ppt_redu, nwFls)]
                if (length(appt_nwdir) > 1) {
                  appt_nwdir <- appt_nwdir[1]
                  print(paste0("Duplicate participant ID:", ppts[p]))
                }

                apptNW <- read.csv(appt_nwdir)
                int <- apptNW$ts[2] - apptNW$ts[1] # Interval between non-wear timestamps

                for (j in 1:nrow(apptNW)){
                  if (apptNW$nonwear[j]){
                    nonWearBool[pptTV > apptNW$ts[j] & pptTV <= apptNW$ts[j]+int] <- TRUE
                  }
                }

                nonWearBool[!pptTVoutslBl] <- FALSE # Make all non-wear inside sleep window = FALSE

                onoffvector_cons[nonWearBool] <- NA # Re-write non-wear times as NA
                onoffvector_sibs[nonWearBool] <- NA
                onoffvector_GGIRnapwind[nonWearBool] <- NA
                onoffvector_GGIRWASOnapwind[nonWearBool] <- NA
                onoffvector_cons_miscl[nonWearBool] <- NA
                onoffvector_sibs_miscl[nonWearBool] <- NA
                onoffvector_GGIRnapwind_miscl[nonWearBool] <- NA
                onoffvector_GGIRWASOnapwind_miscl[nonWearBool] <- NA

              } else {
                print(paste0(ppts[p],
                  ": No non-wear detected. If custom non-wear detection is required (non-GGIR), use 'nonwear_detect' on raw or downsampled accelerometer data to detect non-wear for inclusion in SRI calculation.")
                  )
              }
            }
          }
        }

        # --------------------
        # -> Remove all days with more than 6h NA (overwrite them as NA) ---------
        hto12 <- 12 - originclock
        sec12ts <- originIso + hto12*60*60
        fir12ts <- sec12ts - 24*60*60
        las12ts <- fir12ts + (maxdays + 2)*24*60*60
        setof12s <- as.numeric(seq(from=fir12ts, to=las12ts, by=24*60*60))

        for (j in 1:(length(setof12s)-1)){
          dayBl <- pptTV > setof12s[j] & pptTV <= setof12s[j+1] # Boolean referencing each day

          if (sum(is.na(onoffvector_cons[dayBl])) > 1440/24*exclNAhrs){
            onoffvector_cons[dayBl] <- NA
          }
          if (sum(is.na(onoffvector_sibs[dayBl])) > 1440/24*exclNAhrs){
            onoffvector_sibs[dayBl] <- NA
          }
          if (sum(is.na(onoffvector_GGIRnapwind[dayBl])) > 1440/24*exclNAhrs){
            onoffvector_GGIRnapwind[dayBl] <- NA
          }
          if (sum(is.na(onoffvector_GGIRWASOnapwind[dayBl])) > 1440/24*exclNAhrs){
            onoffvector_GGIRWASOnapwind[dayBl] <- NA
          }
          if (sum(is.na(onoffvector_cons_miscl[dayBl])) > 1440/24*exclNAhrs){
            onoffvector_cons_miscl[dayBl] <- NA
          }
          if (sum(is.na(onoffvector_sibs_miscl[dayBl])) > 1440/24*exclNAhrs){
            onoffvector_sibs_miscl[dayBl] <- NA
          }
          if (sum(is.na(onoffvector_GGIRnapwind_miscl[dayBl])) > 1440/24*exclNAhrs){
            onoffvector_GGIRnapwind_miscl[dayBl] <- NA
          }
          if (sum(is.na(onoffvector_GGIRWASOnapwind_miscl[dayBl])) > 1440/24*exclNAhrs){
            onoffvector_GGIRWASOnapwind_miscl[dayBl] <- NA
          }
        }

        # --------------------
        # -> Calculate SRI: --------
        # Using onset and offset times only ---------
        swv1_1 <- onoffvector_cons[1:(length(onoffvector_cons)-(24*60))]
        swv2_1 <- onoffvector_cons[((24*60)+1):(length(onoffvector_cons))]
        SRI_onoff <- -100 + 200*(1-mean(abs(swv2_1-swv1_1),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

        # Using onset, offset, and SIB gaps >= WASOmin within onset/offset periods ---------
        swv1_2 <- onoffvector_sibs[1:( length(onoffvector_sibs)-(24*60) )]
        swv2_2 <- onoffvector_sibs[((24*60)+1):(length(onoffvector_sibs))]
        SRI_onoffWASO <- -100 + 200*(1-mean(abs(swv2_2-swv1_2),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

        # Using onset, offset + nap window ---------
        swv1_6 <- onoffvector_GGIRnapwind[1:(length(onoffvector_GGIRnapwind)-(24*60))]
        swv2_6 <- onoffvector_GGIRnapwind[((24*60)+1):length(onoffvector_GGIRnapwind)]
        SRI_onoffnap_wind <- -100 + 200*(1-mean(abs(swv2_6-swv1_6),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

        # Using onset, offset + WASO + nap window ---------
        swv1_7 <- onoffvector_GGIRWASOnapwind[1:(length(onoffvector_GGIRWASOnapwind)-(24*60))]
        swv2_7 <- onoffvector_GGIRWASOnapwind[((24*60)+1):length(onoffvector_GGIRWASOnapwind)]
        SRI_onoffWASOnap_wind <- -100 + 200*(1-mean(abs(swv2_7-swv1_7),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

        # As above, but after removal of miscalculated nights ---------
        swv1_8 <- onoffvector_cons_miscl[1:(length(onoffvector_cons_miscl)-(24*60) )]
        swv2_8 <- onoffvector_cons_miscl[((24*60)+1):(length(onoffvector_cons_miscl))]
        misclSRI_onoff <- -100 + 200*(1-mean(abs(swv2_8-swv1_8),na.rm=TRUE))

        swv1_9 <- onoffvector_sibs_miscl[1:( length(onoffvector_sibs_miscl)-(24*60) )]
        swv2_9 <- onoffvector_sibs_miscl[((24*60)+1):(length(onoffvector_sibs_miscl))]
        misclSRI_onoffWASO <- -100 + 200*(1-mean(abs(swv2_9-swv1_9),na.rm=TRUE))

        swv1_10 <- onoffvector_GGIRnapwind_miscl[1:(length(onoffvector_GGIRnapwind_miscl)-(24*60))]
        swv2_10 <- onoffvector_GGIRnapwind_miscl[((24*60)+1):length(onoffvector_GGIRnapwind_miscl)]
        misclSRI_onoffnap_wind <- -100 + 200*(1-mean(abs(swv2_10-swv1_10),na.rm=TRUE))

        swv1_11 <- onoffvector_GGIRWASOnapwind_miscl[1:(length(onoffvector_GGIRWASOnapwind_miscl)-(24*60))]
        swv2_11 <- onoffvector_GGIRWASOnapwind_miscl[((24*60)+1):length(onoffvector_GGIRWASOnapwind_miscl)]
        misclSRI_onoffWASOnap_wind <- -100 + 200*(1-mean(abs(swv2_11-swv1_11),na.rm=TRUE))

        # --------------------
        # -> Check SRIdays for each SRI calculation, find those w/ <= 4 days overlapping data, re-label as NA ---------
        swv1.li <- list(swv1_1, swv1_2, swv1_6, swv1_7, swv1_8, swv1_9, swv1_10, swv1_11) # List all SWVs used for SRI calc
        swv2.li <- list(swv2_1, swv2_2, swv2_6, swv2_7, swv2_8, swv2_9, swv2_10, swv2_11)

        SRIdays.li <- list()
        for (i in 1:length(swv1.li)) {
          SRIdays.li[[i]] <- sum(!is.na(swv1.li[[i]]*swv2.li[[i]]))/24/60
        }

        # --------------------
        # -> Choose required SRI score and 'SRIdays' ----------
        SRIvec <- c(SRI_onoff,SRI_onoffWASO,SRI_onoffnap_wind,SRI_onoffWASOnap_wind,
                    misclSRI_onoff,misclSRI_onoffWASO,misclSRI_onoffnap_wind,misclSRI_onoffWASOnap_wind)

        SRIvec[SRIdays.li < minSRIdays] <- NA # Call all SRI values with fewer than 4 days NA

        t.nap <- c(F,F,T,T,F,F,T,T)
        if (use.naps == FALSE){
          t.nap <- !t.nap
        }
        t.WASO <- c(F,T,F,T,F,T,F,T)

        if (use.WASO == FALSE){
          t.WASO <- !t.WASO
        }

        t.miscal <- c(F,F,F,F,T,T,T,T)
        if (use.miscal == FALSE){
          t.miscal <- !t.miscal
        }
        t.SRI <- t.nap*t.WASO*t.miscal
        t.SRI <- as.logical(t.SRI)
        SRI <- SRIvec[t.SRI]
        SRIdays <- SRIdays.li[[which(t.SRI)]]

        # --------------------
        # -> Calculate and write SWVs for each SRI version ---------
        if (wr.SWV == TRUE){
          # Extract the SWV on/off values as UNIXtime, in the same format as the light data timestamps
          ts1 <- originmin*60 # Get first timestamp

          SWVlist <- list(onoffvector_cons,onoffvector_sibs,onoffvector_GGIRnapwind,onoffvector_GGIRWASOnapwind,
                          onoffvector_cons_miscl,onoffvector_sibs_miscl,onoffvector_GGIRnapwind_miscl,
                          onoffvector_GGIRWASOnapwind_miscl) # Make a list of all SWVs
          names(SWVlist) <- c("onoff", "onoff_WASO", "onoff_nap", "onoff_WASOnap",
                              "mclonoff", "mclonoff_WASO", "mclonoff_nap", "mclonoff_WASOnap")

          sav <- 0
          for (i in 1:length(SWVlist)){
            ts <- sequence(length(SWVlist[[i]]), ts1, 60) # Generate a vector of timestamps

            ts1df <- data.frame(trans = SWVlist[[i]][1],t = ts1) # First timestamp
            tsenddf <- data.frame(trans = SWVlist[[i]][length(SWVlist[[i]])],t = ts[length(ts)]) # Last timestamp

            ont <- ts[which(diff(SWVlist[[i]]) == 1)+1] # Extract on times
            offt <- ts[which(diff(SWVlist[[i]]) == -1)+1] # Extract off times

            if (length(ont) == 0 | length(offt) == 0){
              print(paste0("Error: No sleep-wake transitions are present, ", ppts[p]))
              break
            }

            ondf <- data.frame(trans = 1, t = ont)
            offdf <- data.frame(trans = 0,t = offt)

            NAont <- ts[which(diff(is.na(SWVlist[[i]])) == 1)+1] # Extract NA on times
            if (length(NAont) > 0){
              NAondf <- data.frame(trans = NA,t = NAont)
            } else {
              NAondf <- data.frame(trans = NA,t = NA)
            }

            NAofft <- ts[which(diff(is.na(SWVlist[[i]])) == -1)+1] # Extract NA off times
            if (length(NAofft) > 0){
              NAoffdf <- data.frame(trans = 0,t = NAofft)
              NAoffdf$trans[SWVlist[[i]][which(diff(is.na(SWVlist[[i]])) == -1)+1] == 1] <- 1
            } else {
              NAoffdf <- data.frame(trans = NA,t = NA)
            }

            catdf <- rbind(ts1df,ondf,offdf,NAondf,NAoffdf,tsenddf) # Concatenate
            orcatdf <- catdf[order(catdf$t),] # Order by timestamp
            orcatdf <- orcatdf[!is.na(orcatdf$t),] # Remove all timestamps == NA

            if (orcatdf$t[length(orcatdf$t)] == orcatdf$t[length(orcatdf$t)-1]){ # If duplicate final values, remove last
              orcatdf <- orcatdf[-nrow(orcatdf),]
            }

            if (i == which(t.SRI)){ # Save only the required SWV for raster plot
              s.orcatdf <- orcatdf
              sav <- 1
            }

            writevec <- c(ppts[p], names(SWVlist[i]), as.vector(as.matrix(orcatdf)))
            write.table(t(writevec), SWVfile, sep = ",", col.names = !file.exists(SWVfile), append = TRUE, row.names=FALSE)
          }
        }

        # --------------------
        # -> Read in percentile, compare w/ SRI --------
        SRI_pctl <- quant$X[which(abs(quant$qu - SRI) == min(abs(quant$qu - SRI)))]
        if(length(SRI_pctl) > 1){
          SRI_pctl <- SRI_pctl[1]
        }
        if (is.na(SRI)){
          SRI_pctl <- NA
        }

        # --------------------
        # -> Plot / save raster -----
        if (wr.raster == TRUE){
          if (sav == 1){
            raster_from_SWS(SWS=s.orcatdf,
                            rasdir=rasdir,
                            pptName=ppts[p],
                            tz=tz)
          } else {
            print(paste0("Error: Raster cannot be plotted - no sleep-wake transitions, ", ppts[p]))
          }

        }

        # --------------------
        # -> Write SRI data to file ----------
        writevec <- c(ppts[p],SRI,SRIdays,SRI_pctl,misclnightscount)
        write.table(t(writevec), SRIfile, sep = ",", col.names = !file.exists(SRIfile), append = TRUE, row.names=FALSE)
        print(paste0("SRI extracted: ", ppts[p]))
        # --------------------

      }, error = function(e) {
        print(e)
        error_vec1 <- c(ppts[p],na_vec)
        write.table(t(error_vec1),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                    append = TRUE, row.names=FALSE) # Write Error to SRI.csv
      }, finally = {
        next # Skip to next loop iteration

      }) # Other half of the error-catch function
    } # for ppts close

    }, error = function(e) {
      print(e)
      error_vec1 <- c(studyname,na_vec)
      write.table(t(error_vec1),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                  append = TRUE, row.names=FALSE) # Write Error to SRI.csv
    }, finally = {
      next # Skip to next loop iteration

    }) # Other half of the error-catch function

  } # for k close
  # ----------------------------------------------------
}

# ---------------------------------
# Wrapper -------
#' @title Calculate Sleep Regularity Index from .csv Accelerometer Data
#' @description A wrapper that allows for direct calculation of Sleep Regularity Index (SRI) scores from accelerometer data (.csv format).
#'
#' Runs four functions by default: (a) down-sampling data [ds_accel_csv], (b) extraction of non-wear data [nonwear_detect],
#' (c) predicting sleep-wake timing using GGIR [GGIR_from_csv], and (d) calculating SRI after accounting for
#' fragmented sleep patterns and naps [SRI_from_GGIR].
#'
#' Minimum required inputs: 'acceldir', 'col.timestamp', 'col.accel'.
#'
#' @param acceldir Directory containing raw .csv accelerometer files
#' @param col.timestamp Column of raw .csv files containing timestamp
#' @param col.accel Columns of raw .csv files containing x-y-z accelerometer data e.g., c(1:3)
#' @param alloutdir General output directory, default created if not specified
#' @param dsdir Directory for down-sampled files
#' @param nwdir Directory of non-wear data - specify if using custom (i.e., non-GGIR) non-wear detection. GGIR non-wear data is used preferentially, if available
#' @param outputdir Directory of GGIR output
#' @param rmc.col.time Column of timestamps in down-sampled files
#' @param rmc.col.acc Columns of accelerometer data in down-sampled files
#' @param sdThres Standard deviation threshold for non-wear classification, applied per window
#' @param rngThres Range threshold for non-wear classification, applied per window
#' @param tz Timezone (use OlsonNames() for a list of accepted timezone names
#' @param use.naps Specify whether 'naps' are included in SRI calculation
#' @param use.WASO Specify whether 'wake after sleep onset' (WASO) periods are included in SRI calculation
#' @param use.miscal Specify whether to filter out nights of 'miscalculated' sleep onset/offset timing
#' @param use.GGIRnonwear Specify whether to use GGIR's inbuilt non-wear detection
#' @param use.customnonwear Specify whether to use custom (i.e., non-GGIR) non-wear detection (based on van Hees et. al., 2011)
#' @param nonWearInGGIRsleep Specify whether non-wear periods within GGIR's 'sleep windows' are to be included as non-wear
#' @param wr.SWV Specify whether Sleep-Wake Vectors (SWV) are output to file
#' @param wr.raster Specify whether sleep-wake raster plots are output to file
#' @param minSRIdays Minimum number of days of overlapping data to calculate valid SRI scores
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' SRI_from_accel_csv(acceldir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/106_output/ds_output", col.timestamp = 5, col.accel = c(1:3))
#' SRI_from_accel_csv(acceldir = "~/R/Biobank/GGIR_Version_Testing/105_evensample", col.timestamp = 5, col.accel = c(1:3))
#' SRI_from_accel_csv(acceldir = "~/R/Biobank/Light_Data/127_evenfreq", col.timestamp = 5, col.accel = c(1:3))
SRI_from_accel_csv <- function(acceldir = c(),
                               col.timestamp = c(),
                               col.accel = c(),
                               alloutdir = c(),
                               dsdir = c(),
                               nwdir = c(),
                               outputdir = c(),
                               rmc.col.time = 1,
                               rmc.col.acc = c(2:4),
                               sdThres = 0.12753,
                               rngThres = 0.4905,
                               tz = "UTC",
                               use.naps = TRUE,
                               use.WASO = TRUE,
                               use.miscal = TRUE,
                               use.GGIRnonwear = TRUE,
                               use.customnonwear = TRUE,
                               nonWearInGGIRsleep = FALSE,
                               wr.SWV = TRUE,
                               wr.raster = TRUE,
                               minSRIdays = 5,
                               ...){
  # ---------------------------------
  # Create a separate .csv file containing all specified parameters --------
  # Define all variables input to all subsequent functions # --------
  if (length(acceldir) == 0){
    stop("Error: Specify directory containing accelerometer .csv files")
  }

  if (length(alloutdir) == 0){
    alloutdir <- paste0(acceldir,"_output") # Specify path for all output
    if (!dir.exists(alloutdir)) {
      dir.create(alloutdir) # Create directory if it doesn't exist already
    }
  }

  dsdir <- paste0(alloutdir, "/ds_output") # Directory for down-sampled files
  outputdir <- paste0(alloutdir, "/GGIR_output") # Dir of GGIR output folders
  nwdir <- paste0(alloutdir, "/nw_output")

  # ---------------------------------
  # Run Down-sampling --------
  ds_accel_csv(acceldir = acceldir,
               alloutdir = alloutdir,
               dsdir = dsdir,
               col.timestamp = col.timestamp,
               col.accel = col.accel)

  # ---------------------------------
  # Run Non-wear extraction --------
  if (use.customnonwear){
    nonwear_detect(dsdir = dsdir,
                   nwdir = nwdir,
                   rmc.col.time = rmc.col.time,
                   rmc.col.acc = rmc.col.acc,
                   sdThres = sdThres,
                   rngThres = rngThres)
  }

  # ---------------------------------
  # Run GGIR --------
  GGIR_from_csv(dsdir = dsdir,
                alloutdir = alloutdir,
                outputdir = outputdir,
                rmc.col.acc = rmc.col.acc,
                rmc.col.time = rmc.col.time,
                tz=tz)

  # ---------------------------------
  # Run SRI --------
  SRI_from_GGIR(outputdir = outputdir,
                alloutdir = alloutdir,
                nwdir = nwdir,
                use.naps = use.naps,
                use.WASO = use.WASO,
                use.miscal = use.miscal,
                use.GGIRnonwear = use.GGIRnonwear,
                use.customnonwear = use.customnonwear,
                nonWearInGGIRsleep = nonWearInGGIRsleep,
                wr.SWV = wr.SWV,
                wr.raster = wr.raster,
                minSRIdays = minSRIdays)

  # ----------------------------------------------------
  # Completion message -------
  return("GGIR and SRI analysis complete.")

  # ----------------------------------------------------
}

# ----------------------------------------------------
#' @title Calculate Sleep Regularity Index (SRI) from Binary Sleep-Wake Data
#' @description Calculates SRI from a time series of binary sleep-wake data in .csv format. Column 1 contains values 1=sleep,
#' 0=wake, NA=NA, end=recording end. Column 2 contains UNIX timestamps (origin=1970-01-01) identifying the start of each associated column 1 value.
#' Final timestamp specifies recording end.
#'
#' Minimum required input: binarydir. Specify tz if required (default = "UTC").
#'
#' @param binarydir Directory containing sleep diary data
#' @param tz Timezone (use OlsonNames() for a list of accepted timezone names)
#' @param alloutdir General output directory, default created if not specified
#' @param col.trans Column of sleep-wake transition data
#' @param col.timestamp Column of timestamps
#' @param overwr Specify whether to overwrite previous SRI data
#' @param wr.raster Specify whether to output sleep-wake raster plots
#' @param minSRIdays Minimum number of days of overlapping data to calculate valid SRI scores
#'
#' @return
#' @export
#'
#' @examples
#' SRI_from_binary(binarydir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/diarydata", tz="Portugal")
#' SRI_from_binary(binarydir = "~/R/Biobank/Light_Data/127_evenfreq_output/SWS_output", tz="Australia/Melbourne")
#' SRI_from_binary(binarydir = "~/R/Biobank/Light_Data/127_evenfreq_output/SWS_output_errors")
SRI_from_binary <- function (binarydir = c(),
                             tz = "UTC",
                             alloutdir = c(),
                             col.trans = 1,
                             col.timestamp = 2,
                             overwr = FALSE,
                             wr.raster = TRUE,
                             minSRIdays = 5
                             ){
  # ---------------------------------------
  # Specify directories ---------
  if (length(binarydir) == 0){
    stop("Error: Specify directory containing sleep diary files")
  }
  if (length(alloutdir) == 0){
    alloutdir <- paste0(binarydir,"_output")
  }
  if (!dir.exists(alloutdir)){ # If no output directory, create one
    dir.create(alloutdir)
  }
  if (wr.raster == TRUE){
    rasdir <- paste0(alloutdir,"/raster_output")
    if (!dir.exists(rasdir)){
      dir.create(rasdir) # Create directory
    }
  }

  # Get list of sleep diary files ---------
  file_list <- list.files(binarydir, pattern = "*.csv", full.names = TRUE)
  ppt_list <- list.files(binarydir, pattern = "*.csv")

  # Set up SRI.csv file to write to -----------
  SRIfile <- paste(alloutdir,"SRI.csv",sep="/")
  SRIheader <- c("File", "SRI", "SRI_days", "SRI_pctl")
  if (overwr == FALSE){
    if (!file.exists(SRIfile)) {
      write.table(t(SRIheader),SRIfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
    }
  } else {
    write.table(t(SRIheader),SRIfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
  }
  na_vec <- rep(NA,(length(SRIheader)-1)) # Define vector to fill row for error cases

  # Parameters --------------
  exclNAhrs <- 6 # Hours of NA (i.e., non-wear) data per day above which the entire day will be excluded

  # ---------------------------------------
  # Loop over files, extracting SWVs and SRI -------
  for (k in 1:length(file_list)){
    tryCatch({ # Start of code to catch any error in a loop iteration, write error to SRI file, and skip to next loop iteration
      # --------------------
      # -> Read in data, get filename -------
      studyname <- ppt_list[k] # Names of output folders = filename - directory
      appt.rd <- as.data.frame(data.table::fread(file=file_list[k])) # Read in data for this ppt

      # -> Exclude cases where file format is wrong, write error to file ----------
      if (appt.rd[nrow(appt.rd),2] == appt.rd[nrow(appt.rd)-1,2]){
        appt.rd <- appt.rd[-nrow(appt.rd),]
      }

      if (ncol(appt.rd) != 2 | # If not 2 columns
          length(unique(appt.rd[,1])) > 4 | # If more than 3 unique values in first column
          !all((appt.rd[,2][2:(length(appt.rd[,2]))] - appt.rd[,2][1:(length(appt.rd[,2])-1)]) > 0) # If second column is not sequential
      ){
        error_vec <- c(studyname,na_vec)
        write.table(t(error_vec),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                    append = TRUE, row.names=FALSE) # Write Error to SRI.csv
        print(paste("Incorrect input file format: ",studyname,sep=""))
        next
      }

      # Name d.f. headers
      appt <- data.frame(trans = appt.rd[,col.trans], t = appt.rd[,col.timestamp])

      # -> Round times to the nearest 10 sec ----------
      appt$t.rnd <- round(appt$t,-1)

      # -> Re-construct SWV ---------
      dif <- diff(appt$t.rnd/10) # Count difference between on/off/NA values, 10sec interval for SWV
      SWV <- vector()
      for (i in 1:(nrow(appt)-1)){
        SWV <- c(SWV,rep(appt$trans[i],dif[i]))
      }
      SWV <- as.numeric(SWV)

      # -> Remove all days with more than 6h NA (overwrite them as NA) ---------
      originTs <- appt$t.rnd[1]
      dt <- as.POSIXct(appt$t.rnd[1],origin="1970-01-01",tz=tz)
      originclock <- as.numeric(substr(dt,12,13)) + (as.numeric(substr(dt,15,16)))/60 + (as.numeric(substr(dt,18,19)))/60/60 # 24h clock time of first timestamp
      maxdays <- ceiling((range(appt$t.rnd)[2] - range(appt$t.rnd)[1])/60/60/24) # Max days of data this ppt will have
      pptTV <- seq(from=originTs, by=10, length.out=length(SWV))

      hto12 <- 12 - originclock
      sec12ts <- originTs + hto12*60*60
      fir12ts <- sec12ts - 24*60*60
      las12ts <- fir12ts + (maxdays + 2)*24*60*60
      setof12s <- as.numeric(seq(from=fir12ts, to=las12ts, by=24*60*60))

      for (j in 1:(length(setof12s)-1)){
        dayBl <- pptTV > setof12s[j] & pptTV <= setof12s[j+1] # Boolean referencing each day
        if (sum(is.na(SWV[dayBl])) > 1440/24*exclNAhrs){
          SWV[dayBl] <- NA
        }
      }

      # -> Calculate SRI ---------
      SWV1 <- SWV[1:(length(SWV)-(24*60*6) )] # Remove last 24h
      SWV2 <- SWV[((24*60*6)+1):(length(SWV))] # Remove first 24h
      appt.SRI <- -100 + 200*(1-mean(abs(SWV2-SWV1),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

      # -> Read in percentile, compare w/ SRI --------
      appt.SRI_pctl <- quant$X[which(abs(quant$qu - appt.SRI) == min(abs(quant$qu - appt.SRI)))]

      # -> Number of days used to calculate SRI --------
      appt.SRI_days <- sum(!is.na(SWV1*SWV2))/6/60/24

      if (appt.SRI_days < minSRIdays){
        appt.SRI <- NA
        appt.SRI_pctl <- NA
      }

      # ---------------------------
      # -> Write all data to .csv file -----------
      writevec <- c(studyname, appt.SRI, appt.SRI_days, appt.SRI_pctl)
      write.table(t(writevec), SRIfile, sep = ",", col.names = !file.exists(SRIfile), append = TRUE, row.names=FALSE)
      print(paste0("SRI extracted: ", studyname))

      # ---------------------------
      # -> Plot raster and save to file ----------
      SWS <- appt[,c(1,3)]
      if (wr.raster == TRUE){
        raster_from_SWS(SWS=SWS,
                        rasdir=rasdir,
                        pptName = studyname,
                        tz=tz)
      }

      # ---------------------------
    }, error = function(e) {
      print(e)
      error_vec1 <- c(studyname,na_vec)
      write.table(t(error_vec1),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                  append = TRUE, row.names=FALSE) # Write Error to SRI.csv
    }, finally = {
      next # Skip to next loop iteration
    }) # Error-catch function
  }
  # ---------------------------------------
  # Completion message -------
  return("SRI analysis complete.")

  # ---------------------------------------
}

# ----------------------------------------------------
#' @title Extract Individual SWS Files from single SWV File
#' @description Takes SWV file (single file summary of sleep-wake, generated by 'SRI_from_GGIR')
#' and converts to individual sleep-wake summary (SWS) files. Specify whether output SWS
#' files account for naps, WASO, and miscalculated nights.
#'
#' Minimum required input: SWVfile
#'
#' @param SWVfile Location of sleep-wake vector file
#' @param SWSdir Directory to write individual sleep-wake vector summary (SWS) files
#' @param use.naps Specify whether 'naps' are included in SRI calculation
#' @param use.WASO Specify whether 'wake after sleep onset' (WASO) periods are included in SRI calculation
#' @param use.miscal Specify whether to filter out nights of 'miscalculated' sleep onset/offset timing
#'
#' @return
#' @export
#'
#' @examples
#' SWS_from_SWV(SWVfile="~/R/Biobank/Light_Data/127_evenfreq_output/SWV.csv")
#' SWS_from_SWV(SWVfile="~/R/Biobank/Light_Data/265_evenfreq_output/SWV.csv")
SWS_from_SWV <- function(SWVfile = c(),
                     SWSdir = c(),
                     use.naps = TRUE,
                     use.WASO = TRUE,
                     use.miscal = TRUE
                     ){
  # Check inputs ----------
  if (length(SWVfile) == 0){
    stop("Error: Specify location of sleep-wake vector (SWV) file")
  }
  if (length(SWSdir) == 0){
    basedir <- dirname(SWVfile)
    SWSdir <- paste0(basedir,"/SWS_output")
  }
  if(!dir.exists(SWSdir)){
    dir.create(SWSdir)
  }

  SWVd.coms <- read.csv2(SWVfile)
  SWVd.li <- strsplit(SWVd.coms[[1]],",")

  f2col <- data.frame(ID = rep(NA,length(SWVd.li)), SRItype = NA) # Get first two columns
  for (i in 1:length(SWVd.li)){
    f2col[i,] <- SWVd.li[[i]][1:2]
  }

  # Determine which SRI type should be extracted ----------
  SRIvec <- c("onoff", "onoff_WASO", "onoff_nap", "onoff_WASOnap",
              "mclonoff", "mclonoff_WASO", "mclonoff_nap", "mclonoff_WASOnap")

  t.nap <- c(F,F,T,T,F,F,T,T)
  if (use.naps == FALSE){
    t.nap <- !t.nap
  }
  t.WASO <- c(F,T,F,T,F,T,F,T)
  if (use.WASO == FALSE){
    t.WASO <- !t.WASO
  }
  t.miscal <- c(F,F,F,F,T,T,T,T)
  if (use.miscal == FALSE){
    t.miscal <- !t.miscal
  }
  t.SRI <- t.nap*t.WASO*t.miscal
  t.SRI <- as.logical(t.SRI)

  SRItype <- SRIvec[t.SRI]

  # Extract only SWVs for required SRI type ----------
  SRItype.in <- which(f2col$SRItype == SRItype)
  SWVd.li.SRItype <- SWVd.li[SRItype.in]

  # Re-format as d.f. and write to file ----------
  for (i in 1:length(SWVd.li.SRItype)){

    pptID <- sub(".csv","",SWVd.li.SRItype[[i]][1])
    pptID <- sub(".RData","",pptID) # Extract ID

    pptSWSv <- SWVd.li.SRItype[[i]][-c(1,2)]
    pptSWSdf <- data.frame(trans = pptSWSv[1:(length(pptSWSv)/2)], t = pptSWSv[(length(pptSWSv)/2+1):length(pptSWSv)])
    pptSWSdf$trans[nrow(pptSWSdf)] <- "end"

    pptWrLoc <- paste0(SWSdir, "/", pptID, ".csv")

    add <- 1
    while (file.exists(pptWrLoc)) { # if the file already exists, append to the name
      pptWrLoc <- substr(pptWrLoc, 1, nchar(pptWrLoc)-4)
      pptWrLoc <- paste0(pptWrLoc,"_",add,".csv")
      add <- add+1
    }

    write.table(pptSWSdf, pptWrLoc, col.names=F, row.names = F, sep=",")
  }
}

