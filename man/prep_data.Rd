\name{prep_data} 
\alias{prep_data} 
\usage{prep_data(obesity_file, control_file, ref_name)} 
\title{Loading an preprocessing of the CrescNet dataset.} 
\description{ 
Transform the raw measurements of z-scores. 
} 
\details{The function transform the raw measurements of the initial CrescNet data table to z-scores based on a given reference.} 
\arguments{ 
  \item{obesity_file}{CrescNet file path for obese labeled subjects.} 
  \item{control_file}{CrescNet file path for control labeled subjects.} 
  \item{ref_name}{Which reference to use.} 
} 
\value{ 
 A data.frame contraining all variables of the CrescNet dataset. 
} 
