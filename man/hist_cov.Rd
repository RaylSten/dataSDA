% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hist_cov.R
\name{hist_cov}
\alias{hist_cov}
\title{The covariance of symbolic histogram variable}
\usage{
hist_cov(object, var1, var2, method)
}
\arguments{
\item{object}{A MatH object.}

\item{var1}{A symbolic histogram variable.}

\item{var2}{A symbolic histogram variable.}

\item{method}{The method to calculate covariance. ex: 'BG', 'BD', 'B', 'Wass'}
}
\value{
Return covariance of the histogram variable.
}
\description{
This function compute the covariance of symbolic histogram variable.
}
\examples{
hist_cov(Blood, 'Cholesterol', 'Hemoglobin', method = 'BG')
}
