\name{z_score} 
\alias{z_score} 
\title{Transform raw biometrical measurements to z-scores} 
\usage{z_score(y, L, M, S)} 
\description{ 
Given the LMS parameters and raw measurements for weight, height or BMI,  
the given values are transformed to z-values. 
} 
\details{ 
The function is based on the LMS method of Cole (1991), which is commonly applied for growth curve analysis. 
} 
\arguments{ 
  \item{y}{Vector of values.} 
  \item{L}{Exponent parameter of the Box-Cox function.} 
  \item{M}{The median value.} 
  \item{S}{The standard coefficient of variance.} 
} 
\value{`vector` of transformed values.} 
\references{ 
Cole, TJ, "The LMS method for constructing normalized growth standards", European journal of clinical nutrition 44, 1 (1990), pp. 45-60. 
} 
