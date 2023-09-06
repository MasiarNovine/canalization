\name{transform_data} 
\alias{transform_data} 
\usage{transform_data(cresc_data, ref_name)} 
\title{Loading an preprocessing of the CrescNet dataset.} 
\description{ 
Transform the raw measurements of z-scores. 
} 
\details{The function transform the raw measurements of the initial CrescNet data table to z-scores based on a given reference.} 
\arguments{ 
 \item{cresc_data}{The CrescNet object.} 
 \item{ref_name}{Which reference to use.} 
} 
\value{ 
 A data.frame contraining all variables of the CrescNet dataset. 
} 
