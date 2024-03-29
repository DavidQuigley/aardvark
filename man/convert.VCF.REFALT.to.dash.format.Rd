% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{convert.VCF.REFALT.to.dash.format}
\alias{convert.VCF.REFALT.to.dash.format}
\title{convert VCF-style REF and ALT to equal length with dashes for indels}
\usage{
convert.VCF.REFALT.to.dash.format(s_ref, s_alt)
}
\arguments{
\item{s_ref}{reference nucleotide string}

\item{s_alt}{alternate nucleotide string}
}
\value{
list with values REF, ALT for new formatting
}
\description{
e.g. a REF,ALT pair TTCTT,T -> TTCTT,T----
e.g. a REF,ALT pair T,TTCTT -> T----, TTCTT
e.g. a REF,ALT pair TCT, TAG -> TCT,TAG
}
