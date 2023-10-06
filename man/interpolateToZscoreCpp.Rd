\name{interpolateToZscoreCpp} 
\alias{interpolateToZscoreCpp} 
\title{Linear interpolation of LMS parameters.} 
\usage{interpolateToZscoreCpp(...)} 
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
\seealso{See also \code{\link{getIndicesCpp}}, \code{\link{getProportionCpp}}, \code{\link{interpolateCpp}}.} 
