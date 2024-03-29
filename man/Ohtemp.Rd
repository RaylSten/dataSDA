% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Ohtemp}
\alias{Ohtemp}
\title{30 year trimmed mean daily temperatures interval dataset for the Ohio river basin.}
\format{
A data frame with 161 rows and 7 variables:
\itemize{
    \item \code{ID}: The global historical climatological network (GHCN) station identifier
    \item \code{NAME}: The GHCN station name
    \item \code{STATE}: The two-digit designation for the state in which each station resides
    \item \code{LATITUDE}: Latitude coordinate position
    \item \code{LONGITUDE}: Longitude coordinate position
    \item \code{ELEVATION}: Elevation of the measurement location (meters)
    \item \code{TEMPERATURE}: The 30 year mean daily temperature (tenths of degrees Celsius)
}
}
\source{
\url{https://CRAN.R-project.org/package=intkrige}
}
\usage{
data(Ohtemp)
}
\description{
30 year trimmed mean daily temperatures interval dataset for the Ohio river basin generated from ohtemp dataset.
Intervals are defined by the mean daily maximum and minimum temperatures for
the Ohio river basin from January 1, 1988 - December 31, 2018. The 116
observations in this dataset all had at least 300 daily observations of
temperature in at least 30 of the 31 considered years. The mean was
calculated after trimming 10% of the data in the tails to remove the
influence of potential outliers.
}
\examples{
data(Ohtemp)
}
\keyword{datasets}
