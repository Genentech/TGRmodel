% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GR_curve_fit.R
\name{logistic_4parameters}
\alias{logistic_4parameters}
\title{logistic function for fitting drug-dose response curve}
\usage{
logistic_4parameters(c, Vinf, V0, EC50, h)
}
\arguments{
\item{c}{concentration (can be an array)}

\item{Vinf}{asymptotic value at high concentration}

\item{V0}{asymptotic value at low concentration}

\item{EC50}{mid-point of the curve}

\item{h}{Hill coefficient}
}
\value{
array of response values
}
\description{
\code{logistic_4parameters} returns values based on concentration and fit parameters
}
\details{
returns values based on concentration and fit parameters
}
