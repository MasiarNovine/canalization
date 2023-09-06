\name{get_bmi_differences} 
\alias{get_bmi_differences} 
\title{Calculate the differences of BMI-SDS} 
\description{Cumulative BMI differences of children} 
\details{ 
  For each given age cohort, the differences of consecutive  
  BMI-SDS values of every subject are calculated. If the difference  
  lies in the range given by `age_bounds`, the first of each  
  consecutive value is used to index against the absolute BMI-SDS differences. 
  If `output=list`, the return is a list, with each subject as a list entry.  
  The entries contain the value for the subject ID, the age difference,  
  the BMI-SDS difference and the number of difference values. 
 
  If `output='list-matrix'`, the output is given as a list with the  
  above return but in form of matrices. The values can be plotted,  
  to observe the change of BMI-SDS over the age of individual subjects. 
} 
\arguments{ 
  \item{cresc_data:}{A `data.frame` with CrescNet data.} 
  \item{age_bounds:}{A 2xn `matrix` with age boundaries.} 
  \item{output:}{Either `list`, or `list-matrix`. Default: `list`} 
} 
\value{ 
  Either a list with the subjects and their informations.  
  Otherwise, a list with matrices `matrix`. 
} 
