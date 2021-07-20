# sleepreg
## Calculate Sleep Regularity Index (SRI) scores from accelerometer and/or binary sleep-wake data
### Use cases  
- Use case A: Calculate SRI scores directly from accelerometer data
- Use case B: Calculate SRI scores from GGIR output
- Use case C: Calculate SRI scores from binary sleep-wake data

![Alt text](https://github.com/dpwindred/sleepreg/blob/master/example_flowchart_2.jpg)

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
#### [SRI_from_accel_csv] Calculate Sleep Regularity Index from accelerometer data (.csv format)
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

#### [SRI_from_GGIR] 

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





```
ds_accel_csv(acceldir = c(), 
             alloutdir = c(), 
             dsdir = c(), 
             col.timestamp = c(), 
             col.accel = c())
```

```
GGIR_from_csv(dsdir = c(),
              alloutdir = c(),
              outputdir = c(),
              rmc.col.acc = c(2:4),
              rmc.col.time = 1)
```




```
nonWearDetectDs(dsdir=c(),
                col.timestamp = c(),
                col.accel = c(),
                nwdir=c(),
                sdThres = 0.12753,
                rngThres = 0.4905)
```

```
rollingWindowInd(t=c(),
                 window=c(),
                 step=c())
```


Our package relies upon GGIR (). We acknowledge and thank Vincent van Hees and colleagues for their work in developing this useful package! 
