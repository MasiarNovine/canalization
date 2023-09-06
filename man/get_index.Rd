\name{get_index} 
\alias{get_index} 
\title{Retrieve the pair of indices for linear interpolation.} 
\usage{get_index(age, ref, sex)} 
\description{Retrieve the indices for the interpolation of the LMS parameters} 
\details{ 
The references for the raw measurements only are given for discrete months values, such that  
for ages in-between, the L, M and S parameters must be interpolated. For simplicity, the  
values are linearly interpolated- For this reason, the indices for the lower and upper bound  
in the proximity of the the values of interest must be found by this function. 
} 
\arguments{ 
\item{age}{The age of the child.} 
\item{ref}{The given reference.} 
\item{sex}{The sex of the child.} 
} 
\value{ 
 The relevant indices for interpolation. 
} 
