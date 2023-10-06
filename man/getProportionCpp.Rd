\name{getProportionCpp} 
\alias{getProportionCpp} 
\title{Calculate the proportion factor for precise interpolation.} 
\usage{getProportionCpp(...)} 
\description{ 
 Get the relative proportion, which is necessary for interpolation. 
} 
\details{ 
The references for the raw measurements only are given for discrete 
months values, such that for ages in-between, the L, M and S parameters 
must be interpolated. For simplicity, the values are linearly interpolated. 
For this reason, the indices for the lower and upper bound 
in the proximity of the the values of interest must be found by this 
function. 
} 
\arguments{ 
 \item{...}{Arguments.} 
} 
\value{ 
The proportion of the difference between the lower / upper 
reference age bound and the difference of the the age of the child to the 
upper reference age bound. 
} 
\seealso{See also \code{\link{getIndicesCpp}}} 
