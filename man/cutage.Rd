\name{cutage} 
\alias{cutage} 
\title{Age groups} 
\usage{cutage(x, breaks = c(-Inf, 2, 6, 11, 14, Inf), 
              labels = c("0-2yr", "2-6yr", "6-11yr", "11-14yr", "+14yr"), 
              right = TRUE)} 
\description{Age groups} 
\details{Age groups} 
\arguments{ 
  \item{x}{Age.} 
  \item{breaks}{Age boundaries.} 
  \item{labels}{Age labels.} 
  \item{right}{Whether the intervals are right-closed (default), left-open or both.} 
 } 
\value{ 
  `factor` with age groups. 
 } 
\seealso{\code{\link{cut}}} 
