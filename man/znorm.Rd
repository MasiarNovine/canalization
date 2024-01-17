\name{znorm} 
\alias{znorm} 
\title{Z-normalization} 
\usage{znorm(cresc, ref)} 
\description{ 
Z-transformation of each record 
We need to find the sex-specific indices for each individual record based on the 
interpolated reference table, which for convenience is interpolated in 0.01 years steps. 
} 
\arguments{ 
  \item{cresc}{CrescNet data.table.} 
  \item{ref}{Reference data.table.} 
} 
\value{`data.table` with z-normalized values.} 
