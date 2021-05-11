% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SRI_from_accel_csv_v3b.R
\name{GGIR_from_csv}
\alias{GGIR_from_csv}
\title{Apply GGIR to Down-Sampled Accelerometer Files}
\usage{
GGIR_from_csv(
  dsdir = c(),
  alloutdir = c(),
  outputdir = c(),
  rmc.col.acc = c(2:4),
  rmc.col.time = 1
)
}
\arguments{
\item{dsdir}{Directory for down-sampled files}

\item{alloutdir}{General output directory}

\item{outputdir}{Directory for GGIR output}

\item{rmc.col.acc}{Columns of accelerometer data in down-sampled files}

\item{rmc.col.time}{Column of timestamps in down-sampled files}

\item{rmc.nrow}{Number of accelerometer data rows}
}
\value{

}
\description{
Extracts GGIR output (sleep-wake predictions, sustained inactivity) across all files, specifies GGIR parameters.
}
\examples{
GGIR_from_csv(dsdir = "C:/Users/dan_t/Documents/R/Biobank/SleepRegPackage/106_output/ds_output")
}