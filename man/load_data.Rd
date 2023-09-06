\name{load_data} 
\alias{load_data} 
\title{Loading the CrescNet files.} 
\description{Loading the CrescNet data.} 
\usage{load_data(ob, ctrl)} 
\arguments{ 
  \item{ob}{File path to the CrescNet file containing subjects labeled obese.} 
  \item{ctrl}{File path to the CrescNet file containing subjects labeled control.} 
} 
\details{ 
Based on a `data.table` This function expects a certain file structure given by the two CSV files containing the  
subjects labeled 'obese' and 'control', so initial variable names should not be changed.  
Labeling of the subjects is based on whether any measurements has been over the 95th BMI percentile or not  
using the Kromeyer-Hausschild et al. (2015) reference. The two datasets are pooled together and split into  
a training and a test set. 
} 
\value{A list containing two 'data.table' objects, 'train' and 'test'.} 
\examples{ 
} 
\seealso{ See also \code{\link{build_dataset}} } 
