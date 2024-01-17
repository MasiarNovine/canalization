\name{checkdup} 
\alias{checkdup} 
\title{Check for duplicates} 
\usage{checkdup(cresc)} 
\description{Check for duplicates and keep them, if there 
are harmless otherwise convert to NA to be imputed. 'Harmless' means 
that the values are within 1% of the mean of the duplicates. 
} 
\arguments{ 
  \item{cresc}{CrescNet `data.table`.} 
} 
\value{`data.table` with duplicates removed.} 
