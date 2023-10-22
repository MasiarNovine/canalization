\name{interpolate} 
\alias{interpolate} 
\title{Linear interpolation between an upper and lower value.} 
\description{Perform linear interpolation based on the lower / upper 
reference age bound.} 
\arguments{ 
 \item{prop}{Proportion of the difference value should be used.} 
 \item{low}{Lower value based on the lower reference age bound.} 
 \item{up}{Upper value based on the upper reference age bound.} 
} 
\details{ 
The basic functionality to perform linear interpolation using the 
proportion of the upper and lower age bound, 
relevant for the respective age of a child. 
} 
\value{A 'numeric' value, representing the linear interpolation.} 
\seealso{See also \code{\link{get_indices}}, \code{\link{get_proportion}}, 
\code{\link{interpolate_to_zscore}}} 
