\name{interpolate_reference} 
\alias{interpolate_reference} 
\usage{interpolate_reference(ref_name, years, step, measure, interpol_method, spline_method)} 
\title{Interpolate reference values} 
\description{Takes a reference dataset and interpolates it to a given age range.} 
\details{ 
Based on the reference dataset, the age range and the age step, the function 
interpolates the reference values. Either by linear or spline interpolation. 
} 
\arguments{ 
  \item{ref_name}{A string specifying the reference to be interpolated. 
Possible values are "WHO", "KiGGS" and "Kromeyer". Default: "WHO".} 
  \item{years}{Age range to be interpolated. Default: c(0, 18)} 
  \item{step}{Age step to be used for interpolation. Default: 0.01} 
  \item{measure}{Measure to be interpolated. Either "Weight", "Height" or "BMI".} 
  \item{interpol_method}{Tyope of interpolation. Either "linear" or "spline".} 
  \item{spline_method}{Spline method. See ?spline for available methods. Default: "natural"} 
} 
\value{ Interpolated reference dataset. } 
