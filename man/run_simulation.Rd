\name{run_simulation} 
\alias{run_simulation} 
\title{Run a sampling simulation to emulate CrescNet data filtering} 
\description{ 
Running a simulation based on the general properties of the CrescNet in use. 
} 
\details{ 
  The CrescNet data was filtered by certain selection criteria e.g., the values for 
  the at-risk group have to be over the 97th percentile, while the control 
  contains subjects, which had never any measurement over the 97th percentile. 
  This might introduce a bias into the sample, because those subjects 
  having a high variance of measurement values but do not end obese, 
  would be filtered out of the control, lowering the overall variance of this group. 
  This effect might be small, but must be investigated. 
  The simulation allows for emulate such effects. Note, that currently the age is 
  not modeled. 
} 
\arguments{ 
  \item{cresc_data:}{A `data.frame` with CrescNet data.} 
  \item{age_range:}{A 2 x n `matrix` with age boundaries.} 
} 
\value{ 
  A `list` with as much entries as age cohorts. 
} 
