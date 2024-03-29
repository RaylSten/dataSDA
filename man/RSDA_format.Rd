% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rsda_format.R
\name{RSDA_format}
\alias{RSDA_format}
\title{RSDA Format}
\usage{
RSDA_format(data, sym_type1 = NULL, location = NULL, sym_type2 = NULL, var = NULL)
}
\arguments{
\item{data}{A conventional data.}

\item{sym_type1}{The labels I means an interval variable and $S means set variable.}

\item{location}{The location of the sym_type in the data.}

\item{sym_type2}{The labels I means an interval variable and $S means set variable.}

\item{var}{The name of the symbolic variable in the data.}
}
\value{
Return a dataframe with a label added to the previous column of symbolic variable.
}
\description{
This function changes the format of the data to conform to RSDA format.
}
\examples{
data("Mushroom2")
mushroom.set <- set_variable_format(data = Mushroom2, location = 8, var = "Species")
mushroom.tmp <- RSDA_format(data = mushroom.set, sym_type1 = c("I", "S"),
                            location = c(25, 31), sym_type2 = c("S", "I", "I"),
                            var = c("Species", "Stipe.Length_min", "Stipe.Thickness_min"))
}
