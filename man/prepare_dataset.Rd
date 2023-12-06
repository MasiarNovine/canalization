\name{prepare_dataset} 
\alias{prepare_dataset} 
\usage{prepare_dataset(file_obese, file_control, ref_name, measure, n_train, 
center_age, use_days, encode_sex, use_short_ids, include_baseline, remove_duplicates, seed)} 
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
  \item{file_obese}{CrescNet file path for obese labeled subjects.} 
  \item{file_control}{CrescNet file path for control labeled subjects.} 
  \item{ref_name}{Which reference to use.} 
  \item{measure}{Which measure to use. Either "Weight", "Height" or "BMI".} 
  \item{n_train}{Number of samples used for model building.} 
  \item{center_age}{Should age be centered? Default: FALSE} 
  \item{use_days}{Age given in days, instead of years? Default: FALSE} 
  \item{encode_sex}{Sex encoded as 1 (male) and -1 (female)? Default: FALSE} 
  \item{use_short_ids}{Shorter subject IDs? Default: TRUE} 
  \item{include_baseline}{First measurements as separate variables? Default: TRUE} 
  \item{remove_duplicates}{Remove duplicates? Default: TRUE} 
  \item{seed}{Seed used for random sampling of subjects. Default: 1} 
} 
\value{ 
 A data.frame contraining all variables of the CrescNet dataset. 
} 
