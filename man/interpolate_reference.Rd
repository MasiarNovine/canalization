\name{interpolate_reference} 
\alias{interpolate_reference} 
\usage{interpolate_reference(ref_name, age_range_years, age_step)} 
\title{Interpolate reference values} 
\description{Takes a reference dataset and interpolates it to a given age range.} 
\details{ 
Based on the reference dataset, the age range and the age step, the function 
linearly interpolates the reference values. 
} 
\arguments{ 
  \item{ref_name}{A string specifying the reference to be interpolated. 
Possible values are "WHO", "KiGGS" and "Kromeyer". Default: "WHO".} 
  \item{age_range_years}{Age range to be interpolated. Default: c(0, 18)} 
  \item{age_step}{Age step to be used for interpolation. Default: 0.01} 
} 
\value{ Interpolated reference dataset. } 
