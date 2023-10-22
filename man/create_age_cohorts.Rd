\name{create_age_cohorts} 
\alias{create_age_cohorts} 
\title{Create age cohorts} 
\description{ 
Creating age cohorts. 
} 
\details{ 
  Based on the CrescNet data, different age cohorts are created 
  using given boundaries. The boundaries are given in form of a `matrix`. 
} 
\arguments{ 
  \item{cresc_data}{A `data.frame` with CrescNet data.} 
  \item{age_range}{A 2*n `matrix` with age boundaries.} 
} 
\value{ 
  `list` with as much entries as age cohorts. 
} 
