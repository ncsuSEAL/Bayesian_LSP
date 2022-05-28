% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis_fit.R
\name{PlotAvg}
\alias{PlotAvg}
\title{Visualize the average model fit result. It will show all points as well as the
averaged model fit curve.}
\usage{
PlotAvg(date_vec, vi_vec, avg_fit)
}
\arguments{
\item{date_vec}{the date vector, be sure to convert the vector to "Date"
format or use "yyyy-mm-dd" format string.}

\item{vi_vec}{The vegetation index vector.}

\item{avg_fit}{The model fit object returned by \code{FitAvgModel()}.}
}
\value{
A plot showing the average model fit result.
}
\description{
Visualize the average model fit result. It will show all points as well as the
averaged model fit curve.
}
\examples{
\dontrun{
avg_dt <- FormatAvgData(landsatEVI2$date, landsatEVI2$evi2)
avg_fit <- FitAvgModel(landsatEVI2$date, landsatEVI2$evi2)
PlotAvg(landsatEVI2$date, landsatEVI2$evi2, avg_fit)
}
}