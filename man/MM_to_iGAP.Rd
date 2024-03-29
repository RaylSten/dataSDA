% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MM_to_iGAP.R
\name{MM_to_iGAP}
\alias{MM_to_iGAP}
\title{MM to iGAP}
\usage{
MM_to_iGAP(data)
}
\arguments{
\item{data}{The dataframe with the MM format.}
}
\value{
Return a dataframe with the iGAP format.
}
\description{
To convert MM format to iGAP format.
}
\examples{
data(Face.iGAP)
Face <- iGAP_to_MM(Face.iGAP, 1:6)
MM_to_iGAP(Face)
}
