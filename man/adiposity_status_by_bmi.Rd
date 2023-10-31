\name{adiposity_status_by_bmi} 
\alias{adiposity_status_by_bmi} 
\title{Weight status assignment} 
\description{Assigns a weight status for each observation in the data.} 
\details{ 
Based on given BMI cutoffs, the adiposity status for each observation is determined. 
The argument `cutoffs` must be a named list of 2-element vectors: 
The first element denotes the lower (inclusive) bound, while the second element 
denotes the upper (exclusive) bound. 
} 
\arguments{ 
  \item{cresc_data}{A data.table containing CrescNet data.} 
  \item{ref_name}{Reference name. Either "WHO", "KiGGS" or "Kromeyer-Hauschild". Default: "WHO".} 
  \item{cutoff_bounds}{A named list with BMI cutoffs given as a 2-element vector. Default: list(underweight = c(NA, 18.5), normal = c(18.5, 25), overweight = c(25, 30), obese = c(30, NA))} 
  \item{age_range_years}{Age range in years. Default: c(0, 18)} 
  \item{age_step}{Age step. Default: 0.01} 
} 
\value{A data.table containing an extra column with the assigned weight status for each observation.} 
