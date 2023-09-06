\name{z_score_interpolation} 
\alias{z_score_interpolation} 
\title{Linear interpolation of LMS parameters.} 
\description{ 
 Interpolation of raw values of either weight, height or BMI to z-scores. 
} 
\arguments{ 
 \item{value}{Value to interpolate.} 
 \item{age}{The age of the child.} 
 \item{ref}{The reference data.frame.} 
 \item{type}{The type of the value.} 
} 
\details{ 
For getting z-scores for values for ages, not covered by references tables,  
linear transformation for ages in close proximity are used. 
} 
\value{ 
A 'vector' of z-score transforms. 
} 
\seealso{See also \code{\link{get_index}}, \code{\link{get_proportion}}, \code{\link{lint}}.} 
