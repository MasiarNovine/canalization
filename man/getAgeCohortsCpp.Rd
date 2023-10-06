\name{getAgeCohortsCpp} 
\alias{getAgeCohortsCpp} 
\title{Assign age bounds to observations.} 
\usage{getAgeCohortsCpp(...)} 
\description{ 
Get the age bounds for each observation based on the age of the subject. 
} 
\arguments{ 
 \item{...}{Arguments.} 
} 
\details{ 
For each observation the age is compared to the given age bounds, such that the given age range 
is assigned to the respective observation, if the age at the time of measurement lies in the 
given range. Note, that if the age cannot be assigned to any range, the assigned value is zero. 
} 
\value{ 
A list of with upper and lower age bounds. 
} 
\seealso{See also \code{\link{getIndicesCpp}}, \code{\link{getProportionCpp}}, \code{\link{interpolateCpp}}.} 
