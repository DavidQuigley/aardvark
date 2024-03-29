% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{AW_seq}
\alias{AW_seq}
\title{helper function to return a slice of an AlignmentWindow DNAString between pos_start and pos_end}
\usage{
AW_seq(AW, pos_start, pos_end)
}
\arguments{
\item{AW}{Alignment Window object to examine}

\item{pos_start}{first position to return}

\item{pos_end}{last position to return}
}
\value{
DNAString subsetted by pos_start, pos_end
}
\description{
helper function to return a slice of an AlignmentWindow DNAString between pos_start and pos_end
}
