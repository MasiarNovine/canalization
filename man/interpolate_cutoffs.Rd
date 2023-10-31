\name{interpolate_cutoffs} 
\alias{interpolate_cutoffs} 
\usage{interpolate_cutoffs(ref_name, cutoffs, group_names, age_range_years, age_step)} 
\title{Interpolation of percentile BMI cutoffs} 
\description{Interpolate percentile BMI cutoffs for a given reference.} 
\details{ 
 The function linearly interpolates percentile cutoffs for a given reference based on the 
 International Obesity Task Force (IOTF) cutoffs. The reference is interpolated to the age range 
 defined by `age_range_years` and `age_step`. The function returns a `data.table` with the 
 percentile cutoffs. 
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
\seealso{\code{\link{interpolate_reference}}} 
