% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/READ-read_MEROPS.R
\name{read_merops}
\alias{read_merops}
\title{Read the output of merops and extract abundance profile}
\usage{
read_merops(merops_path, write=FALSE, profile=TRUE)
}
\arguments{
\item{merops_path}{a path where Merops output data are. They
should have the extension .txt and all files in the path are the ones that
need to be read. Output data should have 6 columns with the bin names
followed by the Genes obtained in every algorithm (HMMER,Hotpep,DIAMOND),
column 'Signalp' indcating if a Peptide signal is found and a column
'#ofTools" indicating the number of algorithms that found this Gene.}

\item{write}{a logical value indicating to save the data imported
as a formatted table with .tsv extension with a time stamp and it will be
located in your current workin directory}

\item{profile}{a logical value indicating if you want to print a profile
or not.}
}
\description{
read_merops calculates the abundance of each protease within the
bins based on the merops output from run_merops.
}
\examples{
\dontrun{
read_merops("C:/Users/bins/")
}
}
