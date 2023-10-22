\name{interpolate_cutoff_reference} 
\alias{interpolate_cutoff_reference} 
\usage{interpolate_cutoff_reference(ref_name, cutoffs, group_names, age_range_years, age_step)} 
\title{Interpolation of percentile BMI cutoffs} 
\description{Interpolate percentile BMI cutoffs for a given reference.} 
\details{ 
 The values for the BMI cutoffs are inclusive bounds, i.e. >= 18.5 indicates the start of the normal 
 range and < 25 indicates its end. The group names are , 
 where per default the first group name is "underweight" and the last group name is "obese". 
} 
\arguments{ 
  \item{ref_name}{Reference name. Either "WHO", "KiGGS" or "Kromeyer-Hauschild". Default: "WHO".} 
  \item{cutoffs}{Cutoffs.} 
  \item{group_names}{Names of the weight groups.} 
  \item{age_range_years}{Age range in years.} 
  \item{age_step}{Age step.} 
} 
\value{`data.table` of percentile cutoffs.} 
\references{ 
Cole, MJ, "The LMS method for constructing normalized growth standards", 
European journal of clinical nutrition 44, 1 (1990), pp. 45-60. 
} 
\seealso{ 
\code{\link{interpolate_reference}} 
} 
