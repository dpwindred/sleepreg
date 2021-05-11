
# ----------------------------------------------------
# Takes sleep diary files and calculates:
# - Sleep Regularity Index scores (SRI)
# - Percentile of SRI scores, compared to N ~ 80,000 UK Biobank Participants
# - Mean daily total sleep time (TST)
# - A raster plot of sleep-wake timing

# Notes:
# - Sleep diary files must be in specified format:
#   - Column 1: binary sleep/wake data, 0 = wake start, 1 = sleep start, NA = missing data start
#   - Column 2: times associated with each wake, sleep, or NA start, in seconds.
#   - One file per participant
#   - No header

# D.W. 22/03/21

# ----------------------------------------------------
#' @title Calculate Sleep Regularity Index (SRI) from Sleep Diary Data
#' @description Calculates SRI from sleep-wake data in .csv format, where column 1 contains all sleep-wake transitions (
#' 0 = wake onset, 1 = sleep onset, NA = onset of missing data), column 2 contains UNIX timestamps (seconds, origin = 01/01/1970)
#'
#' @param diarydir Directory containing sleep diary data
#' @param alloutdir General output directory, default created if not specified
#' @param col.trans Column of sleep-wake transition data
#' @param col.timestamp Column of timestamp transition data
#' @param overwr Specify whether to overwrite previous SRI data
#' @param wr.raster Specify whether to output sleep-wake raster plots
#'
#' @return
#' @export
#'
#' @examples
#' SRI_from_diaries(diarydir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/diarydata")
SRI_from_diaries <- function (diarydir = c(), alloutdir = c(), col.trans = 1, col.timestamp = 2, overwr = FALSE, wr.raster = TRUE) {
  # ---------------------------------------
  # Specify directories ---------
  if (length(diarydir) == 0){
    stop("Error: Specify directory containing sleep diary files")
  }
  if (length(alloutdir) == 0){
    alloutdir <- paste0(diarydir,"_output")
  }
  if (!dir.exists(alloutdir)){ # If not output directory, create one
    dir.create(alloutdir)
  }
  if (wr.raster == TRUE){
    rasdir <- paste0(alloutdir,"/raster_output")
    if (!dir.exists(rasdir)){
      dir.create(rasdir) # Create directory
    }
  }

  # Get list of sleep diary files ---------
  file_list <- list.files(diarydir, pattern = "*.csv", full.names = TRUE)

  # Set up SRI.csv file to write to -----------
  SRIfile <- paste(alloutdir,"SRI.csv",sep="/")
  SRIheader <- c("File", "SRI", "SRI_days", "SRI_pctl", "TST")
  if (overwr == FALSE){
    if (!file.exists(SRIfile)) {
      write.table(t(SRIheader),SRIfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
    }
  } else {
    write.table(t(SRIheader),SRIfile,sep=",", col.names=FALSE, row.names=FALSE) #--------------------------------------------------------<<
  }
  na_vec <- rep(NA,(length(SRIheader)-1)) # Define vector to fill row for error cases

  # Record start time -----------
  start_time <- Sys.time()

  # ---------------------------------------
  # Loop over files, extracting SWVs and SRI -------
  for (k in 1:length(file_list)){
    tryCatch({ # Start of code to catch any error in a loop iteration, write error to SRI file, and skip to next loop iteration
      # --------------------
      # -> Read in data, get filename -------
      studyname <- sub(paste(diarydir,"/",sep=""),"",file_list[k])  # Names of output folders = filename - directory
      appt.rd <- as.data.frame(data.table::fread(file=file_list[k])) # Read in data for this ppt

      # -> [FUTURE] Convert timestamps to UNIXsec [outsource in early version] -------
      # -> Exclude cases where file format is wrong, write error to file ----------
      if (ncol(appt.rd) != 2 | # If not 2 columns
          length(unique(appt.rd[,1])) > 3 | # If more than 3 unique values in first column
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

      # -> Determine recording length (days), skip if recording < 5 days ----------
      days <- diff(range(appt$t))/60/60/24
      if (days < 5){ # Don't run this loop if <= 4 nights sleep data
        error_vec <- c(studyname,na_vec)
        write.table(t(error_vec),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                    append = TRUE, row.names=FALSE) # Write Error to SRI.csv
        print(paste("Recording length >=5 days required for accurate SRI calculation: ",studyname,sep=""))
        next
      }

      # -> Round times to the nearest 10 sec ----------
      appt$t.rnd <- round(appt$t,-1)

      # -> Re-construct SWV ---------
      dif <- diff(appt$t.rnd/10) # Count difference between on/off/NA values, 10sec interval for SWV
      SWV <- vector()
      for (i in 1:(nrow(appt)-1)){
        SWV <- c(SWV,rep(appt$trans[i],dif[i]))
      }

      # -> Exclude cases with < 5 nights data after NA exclusion --------
      valsec <- (length(SWV) - sum(is.na(SWV)))*10
      if (valsec < 5*24*60*60){ # Don't run this loop if < 5 nights sleep data
        error_vec <- c(studyname,na_vec)
        write.table(t(error_vec),SRIfile, sep = ",", col.names = !file.exists(SRIfile),
                    append = TRUE, row.names=FALSE) # Write Error to SRI.csv
        print(paste(">=5 days of data required for accurate SRI calculation: ",studyname,sep=""))
        next
      }

      # -> Calculate SRI ---------
      SWV1 <- SWV[1:(length(SWV)-(24*60*6) )] # Remove last 24h
      SWV2 <- SWV[((24*60*6)+1):(length(SWV))] # Remove first 24h
      appt.SRI <- -100 + 200*(1-mean(abs(SWV2-SWV1),na.rm=TRUE)) # Calculate SRI from the two comparison vectors

      # -> Number of days used to calculate SRI --------
      appt.SRI_days <- valsec/60/60/24

      # -> Calculate mean daily TST ---------
      appt.TST <- sum(SWV == 1,na.rm=TRUE)*10/valsec*24

      # ---------------------------
      # -> Plot raster and save to file [grab code from SRI analysis script and turn into a function] ---------
      if (wr.raster == TRUE){
        mdays <- ceiling((range(appt$t)[2] - range(appt$t)[1])/60/60/24) # Max days of data this ppt will have
        appt$tmin <- round(appt$t/60) # Round t to nearest minute

        onind <- which(appt$trans == 1) # Index all sleep onset times
        slt <- vector() # Define empty vector
        for (i in 1:length(onind)){ # Extract all times (minute intervals) of sleep
          slt <- c(slt,appt$tmin[onind[i]]:appt$tmin[(onind[i]+1)])
        }

        sttso <- appt$tmin[2] - appt$tmin[1] # Time difference between recording start and first sleep onset

        rdf <- data.frame(t = slt, tabs = (slt - slt[1] + sttso), day = mdays, tras = (slt - slt[1] + sttso))
        cuts <- seq(from=1440,by=1440,length.out=(mdays-1))
        for (i in 1:length(cuts)){
          rdf$day[rdf$tabs >= cuts[i]] <- rdf$day[rdf$tabs >= cuts[i]] - 1
          rdf$tras[rdf$tabs >= cuts[i]] <- rdf$tras[rdf$tabs >= cuts[i]] - 1440
        }
        rdf$trash <- rdf$tras/60 # In hours

        as.POSIXct(appt$tmin[1]*60,origin="1970-01-01",tz="UTC")
        rot <- as.POSIXct(appt$tmin[1]*60,origin="1970-01-01",tz="UTC")
        oh <- as.numeric(substr(rot,12,13)) + as.numeric(substr(rot,15,16))/60

        rdf$tUTCh <- rdf$trash + oh

        sq <- seq(from=oh,to=(oh+24),by=4)
        sq[sq >= 24] <- sq[sq >= 24] - 24
        abb_x <- sq
        abb_y <- rev(seq(from=max(rdf$day),to=1,by=-1))

        rst <- ggplot2::ggplot(rdf,aes(x=tUTCh, y=day)) +
          geom_point(size=6,shape=15) +
          scale_y_continuous(breaks=seq(from=max(rdf$day),to=1,by=-1),labels=abb_y) +
          scale_x_continuous(breaks=seq(from=oh,to=(oh+24),by=4),limits=c(oh,(oh+24)),labels=abb_x) +
          xlab("Time (UTC, h)") +
          ylab("Day") +
          theme_classic() +
          theme(text = element_text(size = 18), axis.text = element_text(size = 16)) # Change axis height

        w <- 200 # Width in mm
        h <- 7/12*w*(1+5/42*(mdays-7)) # Height (depends on number of days of data)

        ggsave(paste0(rasdir,"/",studyname,".tiff"), device='tiff',
               plot=rst, width=w, height=h, units="mm", dpi=1400,
               limitsize = FALSE)
        print(paste0("Raster plot saved: ", studyname))
      }

      # ---------------------------
      # -> Read in percentile, compare w/ SRI --------
      appt.SRI_pctl <- quan$X[which(abs(quan$qu - appt.SRI) == min(abs(quan$qu - appt.SRI)))]

      # ---------------------------
      # -> Write all data to .csv file -----------
      writevec <- c(studyname, appt.SRI, appt.SRI_days, appt.SRI_pctl, appt.TST)
      write.table(t(writevec), SRIfile, sep = ",", col.names = !file.exists(SRIfile), append = TRUE, row.names=FALSE)
      print(paste0("SRI extracted: ", studyname))

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
  # Completion message and time -------
  end_time1 <- Sys.time()
  t <- end_time1 - start_time
  return(paste("SRI analysis complete. Time = ",t,sep=""))

  # ---------------------------------------
}
# ----------------------------------------------------
# SRI_from_diaries(diarydir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/diarydata",
#                  alloutdir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/outputSRI",
#                  overwr = TRUE)

