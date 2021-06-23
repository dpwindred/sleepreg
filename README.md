# sleepreg
## Calculate Sleep Regularity Index (SRI) scores from accelerometer and/or binary sleep-wake data
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

### Use cases  
This package has been designed for three main use cases:
- Option A: Calculate SRI scores directly from accelerometer data
- Option B: Calculate SRI scores from GGIR output
- Option C: Calculate SRI scores from binary sleep-wake data

Package instructions for each use case:
![Alt text](https://github.com/dpwindred/sleepreg/blob/master/example_flowchart_2.jpg)




Our package relies upon GGIR (). We acknowledge and thank Vincent van Hees and colleagues for their work in developing this exceptionally useful package! 
