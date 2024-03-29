% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aardvark.R
\name{realign_right_softclip}
\alias{realign_right_softclip}
\title{attempt to realign a soft-clipped region on the right side of an alignment}
\usage{
realign_right_softclip(read, min_percent_realigned, align_window)
}
\arguments{
\item{read}{aardvark::Read object on which to attempt the realignment}

\item{min_percent_realigned}{minimum percent of the soft-clipped read segment to accept a new alignment}

\item{align_window}{aardvark::AlignmentWindow object that holds reference sequence search space for realignment}
}
\value{
A modified aardvark::ead object
}
\description{
attempt to realign a soft-clipped region on the right side of an alignment
}
