# sleepreg
## Calculate Sleep Regularity Index (SRI) scores from accelerometer and/or binary sleep-wake data
### Use cases  
- Use case A: Calculate SRI scores directly from accelerometer data
- Use case B: Calculate SRI scores from GGIR output
- Use case C: Calculate SRI scores from binary sleep-wake data

<img src="https://github.com/dpwindred/sleepreg/blob/master/SRI_package_flow7.PNG" width="830" height="650">

### Installation of ‘sleepreg’ package
- Open R
- If you’ve installed GGIR previously, remove any pre-existing version:
```
remove.packages("GGIR")
```
- Close R
- Delete the ‘GGIR’ folder in R/win-library directory
- Re-open R
- Download / install GGIR Version 2.0-0
  - Download ‘GGIR_2.0-0.tar.gz’ from https://cran.r-project.org/src/contrib/Archive/GGIR/
  - Replace [your directory] with directory of the downloaded tar.gz file
  ```
  install.packages(“[your directory]/GGIR_2.0-0.tar.gz", repos = NULL, type="source")
  ```
- Install other dependencies
```
install.packages("ggplot")  
install.packages("data.table")
```
- Install ‘sleepreg’ package’
```
install.packages("devtools")
library(devtools)
auth_token <- "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"
install_github(repo = "dpwindred/sleepreg", auth_token = auth_token)
library(sleepreg)
```
- During installation, R may ask whether you want to update to the newest version of GGIR - say no

### Functions
#### [SRI_from_accel_csv] Calculate Sleep Regularity Index (SRI) from accelerometer data (.csv format)
A wrapper that allows for calculation of Sleep Regularity Index (SRI) scores from accelerometer data (.csv format). Input location of .csv accelerometer files. Function outputs files and folders containing SRI scores, raster plots, sleep-wake vectors, GGIR sleep-wake predictions, miscalculated nights, down-sampled files, and non-wear data.

Runs four functions by default: (a) down-sampling data [ds_accel_csv], (b) extraction of non-wear data [nonwear_detect], (c) predicting sleep-wake timing using GGIR [GGIR_from_csv], and (d) calculating SRI after accounting for fragmented sleep patterns and naps [SRI_from_GGIR].

Minimum required inputs: 'acceldir', 'col.timestamp', 'col.accel'.

```
SRI_from_accel_csv(acceldir = c(),
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
                   ...)
```
|Argument|Description|
|---|---|
| acceldir | Directory containing raw .csv accelerometer files |
|col.timestamp	 |   Column of raw .csv files containing timestamp|
|col.accel	|Columns of raw .csv files containing x-y-z accelerometer data e.g., c(1:3)|
|alloutdir	|General output directory, default created if not specified|
|dsdir	|Directory for down-sampled files|
|nwdir	| Directory of non-wear data - specify if using custom (i.e., non-GGIR) non-wear detection. GGIR non-wear data is used preferentially, if available|
|outputdir	|Directory of GGIR output|
|rmc.col.time	|Column of timestamps in down-sampled files|
|rmc.col.acc	|Columns of accelerometer data in down-sampled files|
|sdThres	|Standard deviation threshold for non-wear classification, applied per window|
|rngThres	|Range threshold for non-wear classification, applied per window|
|tz	|Timezone (use OlsonNames() for a list of accepted timezone names)|
|use.naps	|Specify whether 'naps' are included in SRI calculation|
|use.WASO	|Specify whether 'wake after sleep onset' (WASO) periods are included in SRI calculation|
|use.miscal	|Specify whether to filter out nights of 'miscalculated' sleep onset/offset timing|
|use.GGIRnonwear|	Specify whether to use GGIR's inbuilt non-wear detection|
|use.customnonwear	|Specify whether to use custom (i.e., non-GGIR) non-wear detection (based on van Hees et. al., 2011)|
|nonWearInGGIRsleep	|Specify whether non-wear periods within GGIR's 'sleep windows' are to be included as non-wear|
|wr.SWV	|Specify whether Sleep-Wake Vectors (SWV) are output to file|
|wr.raster	|Specify whether sleep-wake raster plots are output to file|
|minSRIdays	|Minimum number of days of overlapping data to calculate valid SRI scores|

#### [SRI_from_GGIR] Calculate Sleep Regularity Index (SRI) from GGIR Output
Uses sleep windows and sustained inactivity bouts from GGIR output to calculate Sleep Regularity Index scores. Accounts for naps and fragmented sleep by identifying periods of 'wake' during GGIR-defined sleep windows and periods of 'napping' outside GGIR-defined sleep windows. Uses sustained inactivity bouts to exclude days where sleep onset and offset times are likely miscalculated. Runs across all "output_xxx" directories within 'outputdir', accounting for both multi-file and single-file GGIR output structures.

Function outputs are the same as in [SRI_from_accel_csv].

Minimum required inputs: 'outputdir', 'nwdir' if use.customnonwear = TRUE (default)

```
SRI_from_GGIR(outputdir = c(),
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
              minSRIdays = 5)
```
#### [SRI_from_binary] Calculate Sleep Regularity Index (SRI) from binary sleep-wake data
Calculates SRI from a time series of binary sleep-wake summary (SWS) data in .csv format. Column 1 contains values 1=sleep, 0=wake, NA=NA, end=recording end. Column 2 contains UNIX timestamps (origin=1970-01-01) identifying the start of each associated column 1 value. Final timestamp specifies recording end.

Minimum required input: 'binarydir'. Specify 'tz' if required (default = "UTC").
```
SRI_from_binary(binarydir = c(),
                tz = "UTC",
                alloutdir = c(),
                col.trans = 1,
                col.timestamp = 2,
                overwr = FALSE,
                wr.raster = TRUE,
                minSRIdays = 5)
```
|Argument|Description|
|---|---|
|binarydir	|Directory containing sleep diary data|
|tz	|Timezone (use OlsonNames() for a list of accepted timezone names)|
|alloutdir	|General output directory, default created if not specified|
|col.trans	|Column of sleep-wake transition data|
|col.timestamp	|Column of timestamps|
|overwr	|Specify whether to overwrite previous SRI data|
|wr.raster	|Specify whether to output sleep-wake raster plots|
|minSRIdays	|Minimum number of days of overlapping data to calculate valid SRI scores|

#### [raster_from_SWS] Raster plot from sleep-wake summary (SWS) data
Extract raster plots from a binary sleep-wake time series (summarized as a reduced-form data frame). This function is called by [SRI_from_GGIR]. 

Minimum required inputs: 'SWS', 'rasdir', 'pptName'
```
raster_from_SWS(SWS = c(),
                rasdir = c(),
                pptName = c(),
                tz = "UTC")
```
|Argument|Description|
|---|---|
|SWS	|Sleep-wake summary (SWS) data (2 column d.f.)|
|rasdir	|Raster output directory|
|pptName	|Name of participant or file|

#### [ds_accel_csv] Down-sample accelerometer files
Down-samples .csv files to 1Hz, increasing speed and creating required input format for 'GGIR_from_csv' function.

Minimum required inputs: 'acceldir', 'col.timestamp', 'col.accel'.
```
ds_accel_csv(acceldir = c(), 
             alloutdir = c(), 
             dsdir = c(), 
             col.timestamp = c(), 
             col.accel = c())
```
#### [nonwear_detect] Non-wear detection
Evaluates epoch-by-epoch non-wear status of accelerometer devices using method described in van Hees et. al. (2011). 15 minute epochs are evaluated based on surrounding 60 minute windows (centered at 15min), where standard deviation < 13mg and range < 50mg in at least two accelerometer axes is required for non-wear classification.

Minimum required inputs: 'dsdir', 'rmc.col.time', 'rmc.col.acc'
```
nonwear_detect(dsdir=c(),
               nwdir=c(),
               rmc.col.time = 1,
               rmc.col.acc = c(2,3,4),
               sdThres = 0.12753,
               rngThres = 0.4905)
```
#### [GGIR_from_csv] Apply GGIR to down-sampled accelerometer data (.csv format)
Specifies parameters and implements GGIR (Migueles et. al., 2019) across all .csv accelerometer files (frequency = 1Hz) in 'acceldir', extracting sleep-wake predictions and sustained inactivity bouts.

Minimum required input: 'acceldir'
```
GGIR_from_csv(dsdir = c(),
              alloutdir = c(),
              outputdir = c(),
              rmc.col.acc = c(2:4),
              rmc.col.time = 1)
```
#### [SWS_from_SWV] Extract individual SWS files from single SWV file
Takes SWV file (single file summary of sleep-wake, generated by 'SRI_from_GGIR') and converts to individual sleep-wake summary (SWS) files. Specify whether output SWS files account for naps, WASO, and miscalculated nights.

Minimum required input: 'SWVfile'
```
SWS_from_SWV(SWVfile = c(),
             SWSdir = c(),
             use.naps = TRUE,
             use.WASO = TRUE,
             use.miscal = TRUE)
```
|Argument|Description|
|---|---|
|SWVfile	|Location of sleep-wake vector file|
|SWSdir	|Directory to write individual sleep-wake vector summary (SWS) files|

#### [rollingWindowInd] Get rolling window indices
Fits rolling windows to time series and outputs start/end indices of each window.

Note: keep time units consistent between t, window and step. All inputs required. 
```
rollingWindowInd(t=c(),
                 window=c(),
                 step=c())
```
|Argument|Description|
|---|---|
|t	|Time series vector|
|window	|Rolling window size|
|step	|Step size for sliding window along time series|

### Acknowledgements
Our package is built around GGIR (van Hees et. al., 2018; Migueles et. al., 2019). We acknowledge and thank Vincent van Hees and colleagues for their work in developing this useful package! 
