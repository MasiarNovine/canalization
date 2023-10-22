\name{transform_to_centile} 
\alias{transform_to_centile} 
\title{Transformation of z-scores to centiles} 
\usage{transform_to_centile(z, l, m, s)} 
\description{ 
Transform z-scores to centiles 
} 
\arguments{ 
  \item{z}{Vector of z-scores.} 
  \item{l}{Exponent parameter of the Box-Cox function.} 
  \item{m}{The median value.} 
  \item{s}{The standard coefficient of variance.} 
} 
\value{`vector` of transformed values.} 
\references{ 
Cole, MJ, "The LMS method for constructing normalized growth standards", 
European journal of clinical nutrition 44, 1 (1990), pp. 45-60. 
} 
\seealso{ 
\code{\link{transform_to_zscore}} 
} 
