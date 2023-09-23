\name{prep_data} 
\alias{prep_data} 
\usage{prep_data(obesity_file, control_file, ref_name, n_train, 
center_age, age_as_days, encode_sex, seed)} 
\title{Loading an preprocessing of the CrescNet dataset.} 
\description{ 
Transform the raw measurements of z-scores. 
} 
\details{The function transform the raw measurements of the initial 
CrescNet data table to z-scores based on a given reference.} 
\arguments{ 
  \item{obesity_file}{CrescNet file path for obese labeled subjects.} 
  \item{control_file}{CrescNet file path for control labeled subjects.} 
  \item{ref_name}{Which reference to use.} 
  \item{n_train}{Number of samples used for model building.} 
  \item{center_age}{Should age be centered? Default: FALSE} 
  \item{age_as_days}{Age given in days, instead of years? Default: FALSE} 
  \item{encode_sex}{Sex encoded as 1 (male) and -1 (female)? Default: FALSE} 
  \item{seed}{Seed used for random sampling of subjects. Default: 1} 
} 
\value{ 
 A data.frame contraining all variables of the CrescNet dataset. 
} 
