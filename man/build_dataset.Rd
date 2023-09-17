\name{build_dataset} 
\alias{build_dataset} 
\title{Preparation of input data.} 
\description{Preparing the specific CrescNet data from the years 2000-2022.} 
\usage{build_dataset(obesity_file, control_file, n_train, center_age, age_as_days, encode_sex, seed)} 
\arguments{ 
  \item{obesity_file}{CrescNet file path for obese labeled subjects.} 
  \item{control_file}{CrescNet file path for control labeled subjects.} 
  \item{n_train}{Number of samples used for model building.} 
  \item{center_age}{Should age be centered? Default: FALSE} 
  \item{age_as_days}{Should age be given in days as integer, instead of years, i.e. double. Default: FALSE} 
  \item{encode_sex}{Should sex be encoded as zero (male) and one (female)? Default: FALSE} 
  \item{seed}{Value of the seed used for random sampling of subjects. Default: 1} 
} 
\details{ 
The function expects file paths to the two CrescNet data, with the necessary variable names.  
Labeling of the subjects is based on whether any measurement has been over the 95th BMI percentile or not  
using the Kromeyer-Hausschild et al. (2015) reference. The data is split into a training and a test set. 
} 
\value{A list containing two 'data.table' objects, 'train' and 'test'.} 
