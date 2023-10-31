\name{calculateCentileCpp} 
\alias{calculateCentileCpp} 
\title{Single centile calculation.} 
\usage{calculateCentileCpp(...)} 
\description{Calculate a single centile based on given LMS parameters.} 
\details{ 
This function deterines the centile based on the given LMS parameters using 
the LMS method of Cole (1991), which is commonly applied for growth curve analysis 
based on the formula `centile = median * (1 + L * S * z)^(1 / L)`. 
} 
\arguments{ 
 \item{...}{Arguments.} 
} 
\value{A single centile value.} 
\seealso{See also \code{\link{calculateCentilesCpp}}, \code{\link{calculateZScoreCpp}}} 
