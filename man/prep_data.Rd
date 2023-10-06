\name{prep_data} 
\alias{prep_data} 
\usage{prep_data(obesity_file, control_file, ref_name, n_train, 
center_age, age_in_days, encode_sex, short_ids, first_as_variables, seed)} 
\title{Loading an preprocessing of the CrescNet dataset.} 
\description{ 
Transform the raw measurements of z-scores. 
} 
\details{The functions expects file paths to specific CrescNet data with the variables 
'subject_id', 'sex', 'gestational_age', 'age', 'height', 'weight', 'height_sds', 
'bmi', 'bmi_sds' and 'age_group'. Next to raw measurements, Z-scores for 'weight', 
'height' and 'bmi' will be given based on the specified reference. Additionally options include 
the number of subjects used for training and testing the model, next to several options to for 
transforming variables. 
 
} 
\arguments{ 
  \item{obesity_file}{CrescNet file path for obese labeled subjects.} 
  \item{control_file}{CrescNet file path for control labeled subjects.} 
  \item{ref_name}{Which reference to use.} 
  \item{n_train}{Number of samples used for model building.} 
  \item{center_age}{Should age be centered? Default: FALSE} 
  \item{age_in_days}{Age given in days, instead of years? Default: FALSE} 
  \item{encode_sex}{Sex encoded as 1 (male) and -1 (female)? Default: FALSE} 
  \item{short_ids}{Shorter subject IDs? Default: TRUE} 
  \item{first_as_variables}{First measurements as separate variables? Default: TRUE} 
  \item{seed}{Seed used for random sampling of subjects. Default: 1} 
} 
\value{ 
 A data.frame contraining all variables of the CrescNet dataset. 
} 
