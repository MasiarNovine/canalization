\name{z_score} 
\alias{z_score} 
\title{z-score transformation} 
\usage{z_score(y, l, m, s)} 
\description{ 
Transform raw biometrical measurements to z-scores 
} 
\details{ 
Given the LMS parameters and raw measurements for weight, height or BMI,  
the given values are transformed to z-values, based on the  
LMS method of Cole (1991), which is commonly applied for growth curve analysis. 
} 
\arguments{ 
  \item{y}{Vector of values.} 
  \item{l}{Exponent parameter of the Box-Cox function.} 
  \item{m}{The median value.} 
  \item{s}{The standard coefficient of variance.} 
} 
\value{`vector` of transformed values.} 
\references{ 
Cole, TJ, "The LMS method for constructing normalized growth standards", European journal of clinical nutrition 44, 1 (1990), pp. 45-60. 
} 
