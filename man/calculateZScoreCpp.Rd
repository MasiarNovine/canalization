\name{calculateZScoreCpp} 
\alias{calculateZScoreCpp} 
\title{z-score transformation} 
\usage{calculateZScoreCpp(...)} 
\description{Transform raw biometrical measurements to z-scores} 
\details{ 
Given the LMS parameters and raw measurements for weight, height or BMI, 
the given values are transformed to z-values, based on the 
LMS method of Cole (1991), which is commonly applied for 
growth curve analysis. 
} 
\arguments{ 
  \item{...}{Arguments.} 
} 
\value{`vector` of transformed values.} 
\references{Cole, TJ, 
"The LMS method for constructing normalized growth standards", 
European journal of clinical nutrition 44, 1 (1990), pp. 45-60.} 
