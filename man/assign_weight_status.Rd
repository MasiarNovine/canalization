\name{assign_weight_status} 
\alias{assign_weight_status} 
\usage{assign_weight_status(cresc, cutoffs)} 
\title{Weight status assignment} 
\description{Assigns a weight status for each observation in the data.} 
\details{ 
Based on given percentile cutoffs and categories, each observation is categorized. 
The argument `cutoffs` must be a named list of lists: 
The first list denotes the category, while the second list includes following entries 
with the following named entries: `lower`, `upper`, `inclower`, `incupper`, i.e. the 
lower / upper bound (as a numeric value), whether the lower / upper bounds should be 
inclusive or not (as a boolean), respectively. `lower` and `upper` are mandatory, while 
the two other entries are not. If not given, the default assumes inclusive bounds. 
} 
\arguments{ 
  \item{cresc}{A data.table containing CrescNet data.} 
  \item{cutoffs}{BMI cutoffs used. Default: a } 
} 
\value{ 
 A data.table containing an extra column with the assigned weight status for each observation. 
} 
