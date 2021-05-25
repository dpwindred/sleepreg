
# ----------------------------------------------------
# - Calculates SRI based on a) GGIR onset/offset times, b) GGIR onset/offset times and user-defined WASO,
#   c) GGIR onset/offset times and naps based on a rolling window, and d) both (b) and (c)
# - Detects miscalculated GGIR onset/offset times (based on SIBs within a specified pre/post window)
# - Re-calculates all four SRI scores excluding nights with miscalculated times
# --------------------------------------------------
# Notes:
# - Can we add examples of files to the vignettes, or with the package?
# ---------------------------------------------------
# Testing:
# - Converted .csv files from Geneactivs
#   - It's possible to extract as 1hz using the geneactiv software

# - [f] GGIR output from analysis of .bin Geneactiv files (GGIR v2.0-0)

# ---------------------------------------------------
# Additions / fixes:
# - Need to consider timezones / give users an option to specify the timezone in which data was collected. They will give UNIX time,
#   which is absolute relative to UTC, but participants could be anywhere - specify timezone. Only important for the raster plot.
# - Add NAs to raster plots!
# - Make newer versions of GGIR work on the .csv output

# ---------------------------------------------------
# Function structure --------
# - Overall wrapper: SRI_from_accel_csv
#   - Runs: ds_accel_csv (downsamples files)
#   - Runs: GGIR_from_csv (runs GGIR from downsampled files)
#   - Runs: SRI_from_GGIR (extracts SRI from GGIR output)

#   - ds
#     - define fn
#       - define non-def variables
#       - check for / create dir
#       - fn content

#   - ggir
#     - define fn
#       - define non-def variables
#       - install ggir
#       - check for / create dir
#       - fn content

#   - SRI
#     - define fn
#       - define non-def variables
#       - check for / create dir / files
#       - fn content

# - wrapper
#   - define all variables input to all subsequent functions
#   - check for required packages (except GGIR)
#   - run ds fn
#   - run ggir fn
#   - run SRI fn


# ---------------------------------
# Down-sampling --------
#' @title Down-Sample Accelerometer Files
#' @description Down-samples .csv files to 1Hz, increasing speed and creating required input format for 'GGIR_from_csv' function.
#'
#' Minimum required inputs: 'acceldir', 'col.timestamp', 'col.accel'.
#'
#' @param acceldir Directory containing raw .csv accelerometer files
#' @param alloutdir General output directory
#' @param dsdir Directory for down-sampled files
#' @param col.timestamp Column of raw .csv files containing timestamp
#' @param col.accel Columns of raw .csv files containing x-y-z accelerometer data e.g., c(1:3)
#'
#' @return
#' @export
#'
#' @examples
#' ds_accel_csv(acceldir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/106", col.timestamp = 1, col.accel = c(2,3,4))
#'
ds_accel_csv <- function(acceldir = c(), alloutdir = c(), dsdir = c(), col.timestamp = c(), col.accel = c()){
  ### Need to add in error msg here for if columns of timestamp and accelerometer data aren't specified
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
  for (i in 1:length(fl)){ # Loop over files
    wrdir <- paste(dsdir,fl[i],sep="/")
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
# GGIR --------
#' @title Apply GGIR to Down-Sampled Accelerometer Files
#' @description Specifies parameters and implements GGIR (Migueles et. al., 2019) across all .csv accelerometer files (frequency = 1Hz)
#' in 'acceldir', extracting sleep-wake predictions and sustained inactivity bouts.
#'
#' Minimum required inputs: 'acceldir'
#'
#' @param dsdir Directory for down-sampled files
#' @param alloutdir General output directory
#' @param outputdir Directory for GGIR output
#' @param rmc.col.acc Columns of accelerometer data in down-sampled files
#' @param rmc.col.time Column of timestamps in down-sampled files
#' @param rmc.nrow Number of accelerometer data rows
#'
#' @return
#' @export
#'
#' @examples
#'GGIR_from_csv(dsdir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/106_output/ds_output")
GGIR_from_csv <- function(dsdir = c(), alloutdir = c(), outputdir = c(),
                          rmc.col.acc = c(2:4), rmc.col.time = 1){
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
  # # Install GGIR, require ver. 2.0-0, check for and remove any other version --------
  # if(require(GGIR)){
  #   if(sessionInfo()$other$GGIR$Version != "2.0-0" & !is.null(sessionInfo()$other$GGIR$Version)){
  #     uninstall.packages("GGIR") # Remove non-2.0-0 ver. of the pkg
  #     install.packages(GGIRinstdir, repos = NULL, type="source") # Install GGIR version from file
  #   }
  # } else {
  #   install.packages(GGIRinstdir, repos = NULL, type="source") # Install GGIR version from file
  # }
  # library(GGIR)
  # ## Might need to restart the R session from within the loop to make the pkg work??

  # ----------------------------------------------------
  # Run GGIR across downsampled files --------
  file_list <- list.files(dsdir, pattern = "*.csv", full.names = TRUE) # List of downsampled files
  for (k in 1:length(file_list)){
    tryCatch({ # Start of code to catch any error in a loop iteration, write error to SRI file, and skip to next loop iteration
      # --------------------
      # [f] Define 'studyname' -------
      studyname <- sub(paste(dsdir,"/",sep=""),"",file_list[k])  # Names of output folders = filename - directory

      # --------------------
      # [nf] Check if outputdir already exists -------
      checkdir <- paste(outputdir,"/output_",studyname,sep="")
      if (dir.exists(checkdir)) { # If the output folder directory already exists, skip participant
        next
      }
      ## Do we want to not run the script for ppts where the folder already exists -> is there any case where the folder might be partially done?
      ## Can we even overwrite it if that is the case??

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
        rmc.desiredtz = "Europe/London", # Timezone in which device was configured and expriments took place. If experiments took place in a
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
        dofirstpage = FALSE, #First page of the pdf with simple summary histograms
        viewingwindow = 2 #viewing window of the visual report (1 centres at day and 2 centres at night)
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
# SRI --------
#' @title Calculate Sleep Regularity Index (SRI) from GGIR Output
#' @description Uses sleep windows and sustained inactivity bouts from GGIR output to calculate Sleep Regularity Index scores. Accounts for
#' multiphasic and broken sleep by identifying periods of 'wake' during GGIR-defined sleep windows and periods of 'napping' outside
#' GGIR-defined sleep windows. Uses sustained inactivity bouts to exclude days where sleep onset and offset times are likely
#' miscalculated. Runs across all "output_xxx" directories within 'outputdir', accounting for both multi-file and single-file GGIR output
#' structures.
#'
#' Additional outputs: sleep-wake raster plots, sleep variables, summary of miscalculated nights, sleep-wake vectors.
#'
#' Minimum required inputs: 'outputdir'
#'
#' @param outputdir Directory of GGIR output
#' @param alloutdir General output directory
#' @param use.naps Specify whether 'naps' are included in SRI calculation
#' @param use.WASO Specify whether 'wake after sleep onset' (WASO) periods are included in SRI calculation
#' @param use.miscal Specify whether to filter out nights of data likely miscalculated by GGIR
#' @param wr.SWV Specify whether Sleep-Wake Vectors (SWV) are  output to file
#' @param wr.raster Specify whether sleep-wake raster plots are output to file
#'
#' @return
#' @export
#'
#' @examples
#' SRI_from_GGIR(outputdir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/106_output/GGIR_output")
SRI_from_GGIR <- function(outputdir = c(), alloutdir = c(),
                          use.naps = TRUE, use.WASO = TRUE, use.miscal = TRUE, wr.SWV = TRUE, wr.raster = TRUE){
  # ----------------------------------------------------
  # Define undefined variables & check for / create dirs --------
  if (length(outputdir) == 0){
    stop("Error: Specify directory containing GGIR output")
  } # If outputdir not specified, define based on alloutdir
  if (!dir.exists(outputdir)) {
    stop("Error: outputdir doesn't exist")
  } # If no output dir, create one

  if(length(alloutdir) == 0){
    sloca <- unlist(gregexpr("/",outputdir)) # Define alloutdir based on outputdir
    sloc <- sloca[length(sloca)]
    alloutdir <- substr(outputdir,1,(sloc-1))
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
  # Set up file for sleep variables --------
  sleepvarheader <- c("file","night","sleeponset","wakeup","SptDuration","guider_onset","guider_wakeup","guider_SptDuration",
                      "SleepDurationInSpt","calendar_date","SE_GGIR","WASO_GGIR","WASO_Custom","TST_Custom","SE_Custom","naps","miscalculated")
  sleepvarfile <- paste(alloutdir,"sleep_var.csv",sep="/")
  if (!file.exists(sleepvarfile)) {
    write.table(t(sleepvarheader),sleepvarfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
  }

  # ----------------------------------------------------
  # Set up file to write Sleep-Wake Vector (SWV) --------
  SWVheader <- c("file","SRItype")
  SWVfile <- paste(alloutdir,"SWV.csv",sep="/")
  if (!file.exists(SWVfile)) {
    write.table(t(SWVheader),SWVfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
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
  if(length(dir_list) == 0){
    stop("Error: No GGIR output directories available")
  }

  # -----------------------------------------
  na_vec <- rep(NA,(length(SRIheader)-1)) # Define vector to fill row for error cases

  # -----------------------------------------
  # Nap, WASO and miscal parameters -------
  ## Eventually add these to the list of function input params, but need to place restrictions on range
  WASOmin <- 30
  napmin <- 30
  stepsize <- 2 # step size for rolling window (must be a factor of 1440)
  perc <- 0.95 # percentage of window required to be SIB
  misclwindvec <- c(1.5) # windows to check (in hours), applied for all four pre-/post- onset/offset windows. Outputs SRI scores for last choice only.
  misclperc_outsl <- 0.85
  misclperc_insl <- 0.2

  # ----------------------------------------------------
  # Loop over files, extracting SWVs and SRI --------
  for (k in 1:length(dir_list)){
    # a) Assess how many unique participants are contained within this iteration of dir_list. For non-parallel GGIR, will only be one.
    #    For parallel GGIR, will be multiple
    tryCatch({ # Start of code to catch any error in a loop iteration, write error to SRI file, and skip to next loop iteration

      studyname <- sub(paste(outputdir,"/output_",sep=""),"",dir_list[k])  # Names of output folders = filename - directory
      outputfolder <- paste(outputdir,studyname,sep="/output_") # Output folder directory
      nightsumloc <- paste(outputfolder,"/results/QC/part4_nightsummary_sleep_full.csv",sep="") # Nightsummary location
      nightsummary.ppts <- read.csv(nightsumloc, header = TRUE)[,c(1)] # Read ppts in nightsummary
      ppts <- unique(nightsummary.ppts) # List of participants within this GGIR output folder

    for (p in 1:length(ppts)){
      tryCatch({ # Start of code to catch any error in a loop iteration, write error to SRI file, and skip to next loop iteration
        # --------------------
        # Load required data for this ppt -----------
        # Load nightsummary d.f.s
        nightsummary.rd <- read.csv(nightsumloc, header = TRUE)[,c(1,2,3,4,26)] # Read night summary data
        nightsummary <- nightsummary.rd[nightsummary.rd$ID == ppts[p],-1]

        nightsummary2.rd <- read.csv(nightsumloc, header = TRUE)[,c(1,2,3,4,5,7,8,9,14,26)] # Different subset of nightsummary for sl. var. extraction
        nightsummary2 <- nightsummary2.rd[nightsummary2.rd$ID == ppts[p],-1]

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
        # --------------------
        # # [f] Define 'studyname' -------
        # studyname <- sub(paste(outputdir,"/output_",sep=""),"",dir_list[k])  # Names of output folders = filename - directory

        # -----------------------------------------------------------
        # Sleep-wake vector calc
        # --------------------
        # -> [f] Define empty vector, output folder, read nightsummary -----------
        # outputfolder <- paste(outputdir,studyname,sep="/output_") # Define output folder directory
        # nightsumloc <- paste(outputfolder,"/results/QC/part4_nightsummary_sleep_full.csv",sep="") # Define nightsummary location
        # nightsummary <- read.csv(nightsumloc, header = TRUE)[,c(2,3,4,26)] # Read night summary data
        wakeonly_vector <- replicate((nrow(nightsummary)+1)*24*60,0) # Define an empty wake-only vector (all zeros)
        ### Length needs to depend upon length of the file, not specified as just '8'

        # -> [f] Define SWV based on GGIR sleep times only ---------
        nightsummary$onsetmin <- round(nightsummary$sleeponset*60)+(nightsummary$night-1)*24*60 # Convert to min, make continuous
        nightsummary$offsetmin <- round(nightsummary$wakeup*60)+(nightsummary$night-1)*24*60

        nightsummary <- nightsummary[(nightsummary$offsetmin - nightsummary$onsetmin) < 16*60,] # Exclude all nightly 'sleep' times > 16h
        if (nrow(nightsummary) < 1){ # If nightsummary is now empty, skip to next line and write an error.
          error_vec <- c(ppts[p],na_vec,0,"days")
          write.table(t(error_vec),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                      append = TRUE, row.names=FALSE) # Write Error to SRI.csv
          next
        }

        onoffvector_cons <- wakeonly_vector
        for (i in 1:nrow(nightsummary)){ # Make all epochs between onset/offset times equal to one (representing sleep)
          onoffvector_cons[(nightsummary$onsetmin[i]):(nightsummary$offsetmin[i])] <- 1
        }

        # -> [f] Read in sustained inactivity (SIB) data, re-format  --------------------
        # sibloc <- paste(outputfolder,"/meta/ms3.out/",studyname,".RData",sep="")
        # load(sibloc) # Load sustained inactivity data to the GE as "sib.cla.sum"

        sibonoff <- sib.cla.sum[,c("sib.onset.time","sib.end.time")] # Subset sib onset and offset times
        sibonoff <- sibonoff[!sibonoff$sib.onset.time == "",] # Check for and remove all empty sibonoff values

        originmin <- as.numeric(as.POSIXct(nightsummary$calendar_date[1],tz="UTC",tryFormats="%d/%m/%Y"))/60 # Define origin as 00:00 on first listed date

        siboncont <- round((as.numeric(as.POSIXct(sibonoff[,1],tz = "UTC", tryFormats = c("%Y-%m-%dT%H:%M:%S"))))/60)-originmin # Convert each set to UNIXmin
        siboffcont <- round((as.numeric(as.POSIXct(sibonoff[,2],tz = "UTC",tryFormats = c("%Y-%m-%dT%H:%M:%S"))))/60)-originmin
        sibonoff_min1 <- as.data.frame(cbind(siboncont,siboffcont)) # Get a data frame of sib onset/offset in UNIXmin
        sibonoff_min <- sibonoff_min1[sibonoff_min1[,1] > 0,] # Remove any negative value (for SIB onsets outside SRI days)
        sibonoff_min <- sibonoff_min[sibonoff_min[,2] < length(wakeonly_vector),] # Remove all SIB values outside the range of the swv used throughout.

        # -> [f] Define SWV based on GGIR sleep times and WASO (based on SIBs) --------------
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

        # -> [f] Define two SWVs, adding both naps and naps+WASO to GGIR sleep-wake times -----------
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
        # -> Extract sleep variables (nightsummary, custom variables (TST_custom,SE_GGIR,SE_Custom,WASO_GGIR,WASO_Custom,naps)) ---------
        # nightsummary2 <- read.csv(nightsumloc, header = TRUE)[,c(2,3,4,5,7,8,9,14,26)] # Read night summary data
        nightsummary2 <- nightsummary2[nightsummary2$night %in% nightsummary$night,] # Exclude all nights not in 'nightsummary'

        SE_GGIR <- nightsummary2[,"SleepDurationInSpt"] / nightsummary2[,"SptDuration"] # Calculate variables directly based on GGIR output
        WASO_GGIR <- nightsummary2[,"SptDuration"] - nightsummary2[,"SleepDurationInSpt"]

        # Calculate WASO, TST, SE, based on custom WASO calculation method (>WASOmin)
        # -> Bugfix ------
        if (nrow(acceptedsibWASO) >= 1){
          sibmeans <- rowMeans(acceptedsibWASO) # which are greater than offset row n and smaller than onset row n + 1
          fix1acceptedsibWASO <- acceptedsibWASO
          for (i in 1:length(sibmeans)){
            if (sum(sibmeans[i] > nightsummary$offsetmin[-nrow(nightsummary)] & sibmeans[i] < nightsummary$onsetmin[-1]) > 0){
              fix1acceptedsibWASO[i,] <- NA
            }
          }
          FIXEDacceptedsibWASO <- fix1acceptedsibWASO[!is.na(fix1acceptedsibWASO[,1]),]
          # ^^ Had to fix a bug with the acceptedsibWASO data.frame - won't affect anything earlier on in the code though.
        } else {
          FIXEDacceptedsibWASO <- acceptedsibWASO
        }
        # -> Extract WASO and naps ------
        WASOdf <- data.frame(night = nightsummary[,1],WASO = 0) # create empty data.frame for WASO values
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
        # -> Write sleep variable data to file ----------
        sleepvardata <- cbind(rep(ppts[p],nrow(nightsummary2)),nightsummary2,SE_GGIR,WASO_GGIR,WASO_Custom,TST_Custom,SE_Custom,naps,miscalculated)
        write.table(sleepvardata, sleepvarfile, sep = ",", col.names = !file.exists(sleepvarfile), append = TRUE, row.names=FALSE) # Write to sleep variable file

        # --------------------
        # -> [f] Don't calculate SRI for this loop if < 5 nights data -----------
        if (nrow(nightsummary) <= 4){ # Don't calculate SRI for this loop if <= 4 nights sleep data
          error_vec <- c(ppts[p],na_vec)
          write.table(t(error_vec),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                      append = TRUE, row.names=FALSE) # Write Error to SRI.csv
          next
        }

        # --------------------
        # -> [f] Call missing nights 'NA' in SWVs ---------
        missingnights <- setdiff(1:nrow(nightsummary), nightsummary$night)
        if (length(missingnights) > 0){
          for (i in missingnights){
            onoffvector_cons[(721+(i-1)*1440):(720+i*1440)] <- NA
            onoffvector_sibs[(721+(i-1)*1440):(720+i*1440)] <- NA
            onoffvector_GGIRnapwind[(721+(i-1)*1440):(720+i*1440)] <- NA
            onoffvector_GGIRWASOnapwind[(721+(i-1)*1440):(720+i*1440)] <- NA
          }
        }

        # --------------------
        # -> [f] Remove 12h artificial wake time at recording start and end --------------------
        onoffvector_cons_red <- onoffvector_cons[721:(length(onoffvector_cons)-720)] # Define reduced vectors
        onoffvector_sibs_red <- onoffvector_sibs[721:(length(onoffvector_sibs)-720)]
        onoffvector_GGIRnapwind_red <- onoffvector_GGIRnapwind[721:(length(onoffvector_GGIRnapwind)-720)]
        onoffvector_GGIRWASOnapwind_red <- onoffvector_GGIRWASOnapwind[721:(length(onoffvector_GGIRWASOnapwind)-720)]

        # Calculate SRI using onset and offset times only ---------
        swv1_1 <- onoffvector_cons_red[1:(length(onoffvector_cons_red)-(24*60) )] # Remove the first 12h, and the last 36h
        swv2_1 <- onoffvector_cons_red[((24*60)+1):(length(onoffvector_cons_red))] # Remove the first 36h, and the last 12h
        SRI_onoff <- -100 + 200*(1-mean(abs(swv2_1-swv1_1),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

        # Calculate SRI using onset, offset, and SIB gaps >= WASOmin within onset/offset periods ---------
        swv1_2 <- onoffvector_sibs_red[1:( length(onoffvector_sibs_red)-(24*60) )] # Remove the first 12h, and the last 36h
        swv2_2 <- onoffvector_sibs_red[((24*60)+1):(length(onoffvector_sibs_red))] # Remove the first 36h, and the last 12h
        SRI_onoffWASO <- -100 + 200*(1-mean(abs(swv2_2-swv1_2),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

        # Calculate SRI using onset, offset + nap window ---------
        swv1_6 <- onoffvector_GGIRnapwind_red[1:(length(onoffvector_GGIRnapwind_red)-(24*60))] # Remove the first 24h
        swv2_6 <- onoffvector_GGIRnapwind_red[((24*60)+1):length(onoffvector_GGIRnapwind_red)] # Remove the last 24h
        SRI_onoffnap_wind <- -100 + 200*(1-mean(abs(swv2_6-swv1_6),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

        # Calculate SRI using onset, offset + WASO + nap window ---------
        swv1_7 <- onoffvector_GGIRWASOnapwind_red[1:(length(onoffvector_GGIRWASOnapwind_red)-(24*60))] # Remove the first 24h
        swv2_7 <- onoffvector_GGIRWASOnapwind_red[((24*60)+1):length(onoffvector_GGIRWASOnapwind_red)] # Remove the last 24h
        SRI_onoffWASOnap_wind <- -100 + 200*(1-mean(abs(swv2_7-swv1_7),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

        # --------------------
        # -> [f] Number of days used to calculate SRI ----------
        if (length(missingnights) > 0){
          SRIdays <- nrow(nightsummary) - missingnights
        } else {
          SRIdays <- nrow(nightsummary)
        }

        # --------------------
        # -> [f] Re-calculate SRI removing miscalculated onset/offset nights --------
        misclnightscount <- length(misclnights)
        misclSRIdays <- SRIdays - misclnightscount

        onoffvector_cons_miscl <- onoffvector_cons_red
        onoffvector_sibs_miscl <- onoffvector_sibs_red
        onoffvector_GGIRnapwind_miscl <- onoffvector_GGIRnapwind_red
        onoffvector_GGIRWASOnapwind_miscl <- onoffvector_GGIRWASOnapwind_red

        if (length(misclnights) > 0){ # Call all 24h periods with miscalculated sleep onset/offset times 'NA'
          for (i in misclnights){
            onoffvector_cons_miscl[(1+(i-1)*1440):(i*1440)] <- NA
            onoffvector_sibs_miscl[(1+(i-1)*1440):(i*1440)] <- NA
            onoffvector_GGIRnapwind_miscl[(1+(i-1)*1440):(i*1440)] <- NA
            onoffvector_GGIRWASOnapwind_miscl[(1+(i-1)*1440):(i*1440)] <- NA
          }
        }

        if((nrow(nightsummary) - (sum(is.na(onoffvector_cons_miscl))/1440)) >= 5){ # If 5 or more nights of data still exist, calculate SRI
          # Calculate SRI
          swv1_8 <- onoffvector_cons_miscl[1:(length(onoffvector_cons_miscl)-(24*60) )] # Remove the first 12h, and the last 36h
          swv2_8 <- onoffvector_cons_miscl[((24*60)+1):(length(onoffvector_cons_miscl))] # Remove the first 36h, and the last 12h
          misclSRI_onoff <- -100 + 200*(1-mean(abs(swv2_8-swv1_8),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

          swv1_9 <- onoffvector_sibs_miscl[1:( length(onoffvector_sibs_miscl)-(24*60) )] # Remove the first 12h, and the last 36h
          swv2_9 <- onoffvector_sibs_miscl[((24*60)+1):(length(onoffvector_sibs_miscl))] # Remove the first 36h, and the last 12h
          misclSRI_onoffWASO <- -100 + 200*(1-mean(abs(swv2_9-swv1_9),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

          swv1_10 <- onoffvector_GGIRnapwind_miscl[1:(length(onoffvector_GGIRnapwind_miscl)-(24*60))] # Remove the first 24h
          swv2_10 <- onoffvector_GGIRnapwind_miscl[((24*60)+1):length(onoffvector_GGIRnapwind_miscl)] # Remove the last 24h
          misclSRI_onoffnap_wind <- -100 + 200*(1-mean(abs(swv2_10-swv1_10),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

          swv1_11 <- onoffvector_GGIRWASOnapwind_miscl[1:(length(onoffvector_GGIRWASOnapwind_miscl)-(24*60))] # Remove the first 24h
          swv2_11 <- onoffvector_GGIRWASOnapwind_miscl[((24*60)+1):length(onoffvector_GGIRWASOnapwind_miscl)] # Remove the last 24h
          misclSRI_onoffWASOnap_wind <- -100 + 200*(1-mean(abs(swv2_11-swv1_11),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

        } else {
          misclSRI_onoff <- NA
          misclSRI_onoffWASO <- NA
          misclSRI_onoffnap_wind <- NA
          misclSRI_onoffWASOnap_wind <- NA
        }

        # --------------------
        # -> Choose required SRI score ----------
        SRIvec <- c(SRI_onoff,SRI_onoffWASO,SRI_onoffnap_wind,SRI_onoffWASOnap_wind,
                    misclSRI_onoff,misclSRI_onoffWASO,misclSRI_onoffnap_wind,misclSRI_onoffWASOnap_wind)
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
        } else {
          SRIdays <- misclSRIdays
        }
        t.SRI <- t.nap*t.WASO*t.miscal
        t.SRI <- as.logical(t.SRI)
        SRI <- SRIvec[t.SRI]

        # --------------------
        # -> [f] Calculate SWVs for each SRI version ---------
        if (wr.SWV == TRUE){
          # Extract the SWV on/off values as UNIXtime, in the same format as the light data timestamps
          ts1 <- originmin*60 + 60*60*12 # Get first timestamp
          SWVlist <- list(onoffvector_cons_red,onoffvector_sibs_red,onoffvector_GGIRnapwind_red,onoffvector_GGIRWASOnapwind_red,
                          onoffvector_cons_miscl,onoffvector_sibs_miscl,onoffvector_GGIRnapwind_miscl,
                          onoffvector_GGIRWASOnapwind_miscl) # Make a list of all SWVs
          names(SWVlist) <- c("onoff", "onoff_WASO", "onoff_nap", "onoff_WASOnap",
                              "mclonoff", "mclonoff_WASO", "mclonoff_nap", "mclonoff_WASOnap")

          for (i in 1:length(SWVlist)){
            ts <- sequence(length(SWVlist[[i]]), ts1, 60) # Generate a vector of timestamps

            ts1df <- data.frame(trans = SWVlist[[i]][1],t = ts1) # First timestamp
            tsenddf <- data.frame(trans = SWVlist[[i]][length(SWVlist[[i]])],t = ts[length(ts)]) # Last timestamp

            ont <- ts[which(diff(SWVlist[[i]]) == 1)+1] # Extract on times
            ondf <- data.frame(trans = 1,t = ont)

            offt <- ts[which(diff(SWVlist[[i]]) == -1)+1] # Extract off times
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

            if (i == which(t.SRI)){ # Save only the required SWV for raster plot
              s.orcatdf <- orcatdf
            }

            writevec <- c(ppts[p], names(SWVlist[i]), as.vector(as.matrix(orcatdf)))
            write.table(t(writevec), SWVfile, sep = ",", col.names = !file.exists(SWVfile), append = TRUE, row.names=FALSE)
          }
        }

        # --------------------
        # -> Read in percentile, compare w/ SRI --------
        # quan <- read.csv(quantdir) # Read in quantile data
        SRI_pctl <- quan$X[which(abs(quan$qu - SRI) == min(abs(quan$qu - SRI)))]
        if (is.na(SRI)){
          SRI_pctl <- NA
        }

        # --------------------
        # -> Plot /save raster -----
        if (wr.raster == TRUE){
          mdays <- ceiling((range(s.orcatdf$t)[2] - range(s.orcatdf$t)[1])/60/60/24) # Max days of data this ppt will have
          s.orcatdf$tmin <- round(s.orcatdf$t/60) # Round t to nearest minute

          onind <- which(s.orcatdf$trans == 1) # Index all sleep onset times
          slt <- vector() # Define empty vector
          for (i in 1:length(onind)){ # Extract all times (minute intervals) of sleep
            slt <- c(slt,s.orcatdf$tmin[onind[i]]:s.orcatdf$tmin[(onind[i]+1)])
          }

          sttso <- s.orcatdf$tmin[2] - s.orcatdf$tmin[1] # Time difference between recording start and first sleep onset

          rdf <- data.frame(t = slt, tabs = (slt - slt[1] + sttso), day = mdays, tras = (slt - slt[1] + sttso))
          cuts <- seq(from=1440,by=1440,length.out=(mdays-1))
          for (i in 1:length(cuts)){
            rdf$day[rdf$tabs >= cuts[i]] <- rdf$day[rdf$tabs >= cuts[i]] - 1
            rdf$tras[rdf$tabs >= cuts[i]] <- rdf$tras[rdf$tabs >= cuts[i]] - 1440
          }
          rdf$trash <- rdf$tras/60 # In hours

          as.POSIXct(s.orcatdf$tmin[1]*60,origin="1970-01-01",tz="UTC")
          rot <- as.POSIXct(s.orcatdf$tmin[1]*60,origin="1970-01-01",tz="UTC")
          oh <- as.numeric(substr(rot,12,13)) + as.numeric(substr(rot,15,16))/60

          rdf$tUTCh <- rdf$trash + oh

          sq <- seq(from=oh,to=(oh+24),by=4)
          sq[sq >= 24] <- sq[sq >= 24] - 24
          abb_x <- sq
          abb_y <- rev(seq(from=max(rdf$day),to=1,by=-1))

          rst <- ggplot(rdf,aes(x=tUTCh, y=day)) +
            geom_point(size=6,shape=15) +
            scale_y_continuous(breaks=seq(from=max(rdf$day),to=1,by=-1),labels=abb_y) +
            scale_x_continuous(breaks=seq(from=oh,to=(oh+24),by=4),limits=c(oh,(oh+24)),labels=abb_x) +
            xlab("Time (UTC, h)") +
            ylab("Day") +
            theme_classic() +
            theme(text = element_text(size = 18), axis.text = element_text(size = 16)) # Change axis height

          w <- 200 # Width in mm
          h <- 7/12*w*(1+5/42*(mdays-7)) # Height (depends on number of days of data)

          ggsave(paste0(rasdir,"/",ppts[p],".tiff"), device='tiff',
                 plot=rst, width=w, height=h, units="mm", dpi=1400,
                 limitsize = FALSE)
          print(paste0("Raster plot saved: ", ppts[p]))
        }

        # --------------------
        # -> Write SRI data to file ----------
        writevec <- c(ppts[p],SRI,SRIdays,SRI_pctl,misclnightscount)
        write.table(t(writevec), SRIfile, sep = ",", col.names = !file.exists(SRIfile), append = TRUE, row.names=FALSE)
        print(paste0("SRI extracted: ", ppts[p]))

      }, error = function(e) {
        print(e)
        error_vec1 <- c(ppts[p],na_vec,"Error")
        write.table(t(error_vec1),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                    append = TRUE, row.names=FALSE) # Write Error to SRI.csv
      }, finally = {
        next # Skip to next loop iteration

      }) # Other half of the error-catch function
    } # for ppts close

    }, error = function(e) {
      print(e)
      error_vec1 <- c(studyname,na_vec,"Error")
      write.table(t(error_vec1),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                  append = TRUE, row.names=FALSE) # Write Error to SRI.csv
    }, finally = {
      next # Skip to next loop iteration

    }) # Other half of the error-catch function

  } # for k close
  # ----------------------------------------------------
} # Function close

# ---------------------------------
# Wrapper -------
#' @title Calculate Sleep Regularity Index from .csv Accelerometer Data
#' @description A wrapper that allows for direct calculation of Sleep Regularity Index (SRI) scores from accelerometer data (.csv format).
#'
#' Runs three functions: (a) down-sampling data [ds_accel_csv], (b) predicting sleep-wake timing using GGIR [GGIR_from_csv],
#' and (c) calculating SRI after accounting for broken sleep patterns [SRI_from_GGIR].
#'
#' Minimum required inputs: 'acceldir', 'col.timestamp', 'col.accel'.
#'
#' @param acceldir Directory containing raw .csv accelerometer files
#' @param alloutdir General output directory
#' @param dsdir Directory for down-sampled accelerometer files
#' @param col.timestamp Column of raw .csv files containing timestamp
#' @param col.accel Columns of raw .csv files containing x-y-z accelerometer data e.g., c(1:3)
#' @param outputdir Directory for GGIR output
#' @param rmc.col.acc Columns of accelerometer data in down-sampled files
#' @param rmc.col.time Column of timestamps in down-sampled files
#' @param rmc.nrow Number of accelerometer data rows
#' @param use.naps Specify whether 'naps' are included in SRI calculation
#' @param use.WASO Specify whether 'wake after sleep onset' (WASO) periods are included in SRI calculation
#' @param use.miscal Specify whether to filter out nights of 'miscalculated' sleep onset/offset timing
#' @param wr.SWV Specify whether Sleep-Wake Vectors (SWV) are output to file
#' @param wr.raster Specify whether sleep-wake raster plots are output to file
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' SRI_from_accel_csv(acceldir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/106_output/ds_output", col.timestamp = 5, col.accel = c(1:3))
SRI_from_accel_csv <- function(acceldir = c(),
                               alloutdir = c(),
                               dsdir = c(),
                               col.timestamp = c(),
                               col.accel = c(),
                               outputdir = c(),
                               rmc.col.acc = c(2:4),
                               rmc.col.time = 1,
                               rmc.nrow = c(),
                               use.naps = TRUE,
                               use.WASO = TRUE,
                               use.miscal = TRUE,
                               wr.SWV = TRUE,
                               wr.raster = TRUE,
                               ...){
  # ---------------------------------
  # Record start time -----------
  start_time <- Sys.time()

  # Create a separate .csv file containing all specified parameters --------
  # Define all variables input to all subsequent functions # --------
  if (length(acceldir) == 0){
    stop("Error: Specify directory containing accelerometer .csv files")
  }

  alloutdir <- paste0(acceldir,"_output") # Specify path for all output
  if (!dir.exists(alloutdir)) {
    dir.create(alloutdir) # Create directory if it doesn't exist already
  }

  dsdir <- paste0(alloutdir, "/ds_output") # Directory for down-sampled files
  outputdir <- paste0(alloutdir, "/GGIR_output") # Dir of GGIR output folders

  # quantdir <- "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/quantiles_1000.csv" # Location of SRI quantile data
  # GGIRinstdir <- "C:/Users/dan_t/Documents/R/Biobank/GGIR_Version_Testing/Versions/GGIR_2.0-0.tar.gz" # Location of GGIR install files


  # # Check for and load/install required packages (except GGIR) # --------
  # if(!require(readr)){
  #   install.packages("readr")
  # }
  # library(readr) # Or maybe we need to make this data.table / fread??
  #
  # if(!require(data.table)){
  #   install.packages("data.table")
  # }
  # library(data.table)
  #
  # if(!require(installr)){
  #   install.packages("installr")
  # }
  # library(installr)
  #
  # if(!require(ggplot2)){
  #   install.packages("ggplot2")
  # }
  # library(ggplot2)


  # ---------------------------------
  # Run Down-sampling --------
  ds_accel_csv(acceldir=acceldir, alloutdir=alloutdir, dsdir=dsdir, col.timestamp = col.timestamp, col.accel = col.accel)

  # ---------------------------------
  # Run GGIR --------
  GGIR_from_csv(dsdir = dsdir, alloutdir = alloutdir, outputdir = outputdir)

  # ---------------------------------
  # Run SRI --------
  SRI_from_GGIR(alloutdir=alloutdir, outputdir=outputdir)

  # ----------------------------------------------------
  # [nv] Completion message and time -------
  end_time1 <- Sys.time()
  t <- end_time1 - start_time
  return(paste("GGIR and SRI analysis complete. Time = ",t,sep=""))

  # ----------------------------------------------------
}
# ----------------------------------------------------

