\name{interpolate_to_z_score_vector_cpp} 
\alias{interpolate_to_z_score_vector_cpp} 
\title{Linear interpolation of LMS parameters.} 
\usage{interpolate_to_z_score_vector_cpp(...)} 
\description{ 
 Interpolation of a vector of raw values of either weight, height or BMI to z-scores. 
} 
\arguments{ 
 \item{...}{Arguments.} 
} 
\details{ 
For getting z-scores for values for ages, not covered by references tables, 
linear transformation for ages in close proximity are used. 
} 
\value{ 
A 'vector' of z-score transforms. 
} 
\seealso{See also \code{\link{get_indices_cpp}}, \code{\link{get_proportion_cpp}}, \code{\link{lint_cpp}}.} 
