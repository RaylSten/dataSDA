% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iGAP_to_MM.R
\name{iGAP_to_MM}
\alias{iGAP_to_MM}
\title{iGAP to MM}
\usage{
iGAP_to_MM(data, location)
}
\arguments{
\item{data}{The dataframe with the iGAP format.}

\item{location}{The location of the symbolic variable in the data.}
}
\value{
Return a dataframe with the MM format.
}
\description{
To convert iGAP format to MM format.
}
\examples{
data(Abalone.iGAP)
Abalone <- iGAP_to_MM(Abalone.iGAP, 1:7)
}
