% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CALC-calc_percentage.R
\name{calc_percentage}
\alias{calc_percentage}
\title{Calculate the percentage.}
\usage{
calc_percentage (tibble_ko, y_axis, data_experiment = NULL)
}
\arguments{
\item{tibble_ko}{a tibble object from mapping_ko.}

\item{y_axis}{a character, indicating the pathway to analyze.}

\item{data_experiment}{optional. a data frame object containing metadata
information.}
}
\value{
A data frame with the calculated percentage of genes in a certain
pathway of the data.
}
\description{
Calculate the percentage of a certain pathway in the data.
}
\details{
Calculate the percentage of a certain pathway in the data.
}
\examples{
# calc_percentage (ko_bin_mapp, Pathway)   
}
