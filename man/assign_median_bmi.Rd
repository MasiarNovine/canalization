\name{assign_median_bmi} 
\alias{assign_median_bmi} 
\usage{assign_median_bmi(cresc_data, ref_name)} 
\title{Assign median BMI} 
\description{Assign median BMI for a given reference.} 
\details{ 
 Assigning the interpolated BMI median for a given reference for 
 age and sex adjustment. 
} 
\arguments{ 
  \item{cresc_data}{Input data.frame.} 
  \item{ref_name}{Reference name. Either "WHO", "KiGGS" or "Kromeyer-Hauschild". Default: "WHO".} 
} 
\value{`data.table` with BMI median.} 
