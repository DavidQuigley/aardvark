% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aardvark.R
\name{realign_read}
\alias{realign_read}
\title{attempt to improve the alignment of a read through local realignment}
\usage{
realign_read(
  read,
  align_window,
  pathogenic_mutation,
  gr_pathogenic = NULL,
  allow_insertions_in_realign = TRUE,
  min_nt_for_distant_realign = 15,
  min_percent_realigned = 0.9,
  near_bound = 6,
  min_nt_qual = 20
)
}
\arguments{
\item{read}{aardvark::Read object under consideration}

\item{align_window}{aardvark::AlignmentWindow object that holds reference sequence search space for realignment}

\item{pathogenic_mutation}{the pathogenic mutation}

\item{gr_pathogenic}{GenomicRanges spanning the pathogenic mutation}

\item{allow_insertions_in_realign}{whether to allow insertions in the realignment}

\item{min_nt_for_distant_realign}{minimum size for a local realignment that moves a soft-clipped segment to a distant locus}

\item{min_percent_realigned}{minimum percentage of the total length of a soft-clipped segment that must be aligned perfectly}

\item{near_bound}{number of nucleotides adjacent to a germline mutation to consider off of the end of a read for local realignment}

\item{min_nt_qual}{minimal nucleotide quality score to count when evaluating perfect matches in a candidate alignment}
}
\value{
Returns the updated read object
}
\description{
Pass both gr_pathogenic and pathogenic_mutation because making a GenomicRanges
is expensive if you're doing it thousands of times.
}
