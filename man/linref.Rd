\name{linref} 
\alias{linref} 
\title{Linear interpolation of LMS parameters.} 
\description{ 
 Interpolation of raw values of either weight, height or BMI to z-scores. 
 
The references for the raw measurements only are given for discrete 
months values, such that for ages in-between, the L, M and S parameters 
are linearly interpolated. 
 
For getting z-scores for values for ages, not covered by references tables, 
linear transformation for ages in close proximity are used. 
} 
\arguments{ 
 \item{ref}{The reference data.frame.} 
} 
