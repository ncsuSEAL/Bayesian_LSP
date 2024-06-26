% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_fit.R
\name{BLSPFitted}
\alias{BLSPFitted}
\title{Get the BLSP fitted curves}
\usage{
BLSPFitted(blsp_fit, asCI = FALSE)
}
\arguments{
\item{blsp_fit}{The "BlspFit" object.}

\item{asCI}{Logical. Default is \code{TRUE}, which means the fitted curves will be
summarized as median and 95\% credibel interval.}
}
\value{
The fitted curves, or the summary, for each year.
}
\description{
The funciton uses the "BlspFit" object returned by the "FitBLSP" function to
predict the fitted daily VI curves.
}
