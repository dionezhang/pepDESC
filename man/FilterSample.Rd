% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FilterSample.R
\name{FilterSample}
\alias{FilterSample}
\title{Title FilterSample
This function removes outlier samples according to the identification and total abundance.}
\usage{
FilterSample(Peptide)
}
\arguments{
\item{Peptide}{Quantitative peptide processed without untrusted data or contaminant data}
}
\value{
A data list of peptide information without outlier samples
}
\description{
Title FilterSample
This function removes outlier samples according to the identification and total abundance.
}
\examples{
Peptide<-Initiation(test_raw)
Peptide<-FilterData(Peptide,type=test_type,'A','B')
Peptide<-FilterSample(Peptide)
}
