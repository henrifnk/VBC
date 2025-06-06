% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_data.R
\name{map_delta}
\alias{map_delta}
\title{Delta mapping}
\usage{
map_delta(mp, mph, mp_kde, mc_kde, xmin = NA)
}
\arguments{
\item{mp}{\link{double}\cr
uncorrected climate variable in projection period.}

\item{mph}{\link{double}\cr
corrected climate variable in projection period.}

\item{mp_kde}{\link[kde1d:kde1d]{kde1d::kde1d}\cr
a kernel density estimation of the climate variable in projection period.}

\item{mc_kde}{\link[kde1d:kde1d]{kde1d::kde1d}\cr
a kernel density estimation of the climate variable in calibration period.}

\item{xmin}{double(1)\cr
A vector indicating if xmin is a ratio variable type or any other type (NA).}
}
\value{
A climate variable that is corrected by the climate trend in the
model between correction and projection period.
}
\description{
Map the rank delta in model between calibration period \code{mc} and projection
period \code{mp} on corrected projection period \code{mph} to correct for climate
trends. For climate variables bounded at 0, a multiplicative delta scheme is
applied when the projection quantiles are lower than the calibration
quantiles.
}
