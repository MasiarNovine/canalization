\name{capfactorDT} 
\alias{capfactorDT} 
\title{Capitalize the first letter of each word, but keeping factors as factors} 
\usage{capfactorDT(dt, colnames)} 
\description{Capitalize the first letter of each word, but keeping factors as factors} 
\arguments{ 
  \item{dt}{The data.table containing the columns.} 
  \item{colnames}{A character vector with the names for factor columns.} 
} 
\details{ 
This function uses the \code{capwords} function and capitalizes the first letter of 
each word in the given columns of a data.table based on the given column names, 
but will keep factors as factors. 
} 
\value{The initial data.table with the updated column} 
\seealso{See also \code{\link{capwords}}.} 
