\name{ggvario} 
\alias{ggvario} 
\title{Plot a variogram} 
\description{Uses ggplot to plot a variogram} 
\usage{ggvario(v)} 
\arguments{ 
  \item{v}{A data.frame object with columns `dist` and `variog`.} 
} 
\details{ 
The column `dist` contains the distances, while `variog` is the the result 
of the variogram function. 
} 
\value{A `ggplot` object.} 
