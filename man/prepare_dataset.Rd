\name{prepare_dataset} 
\alias{prepare_dataset} 
\usage{prepare_dataset(file_adiposity, file_control, n_train, 
center_age, use_days, encode_sex, short_ids, first_extra, check_dups, 
seed, age_range)} 
\title{Loading and preprocessing of the CrescNet dataset.} 
\description{Prepare the CrescNet dataset.} 
\details{The functions expects file paths to specific CrescNet data with the variables 
'subject_id', 'sex', 'gestational_age', 'age', 'height', 'weight', 'height_sds', 
'bmi', 'bmi_sds' and 'age_group'. `check_dups` will look for duplicated values, 
if their values are indistinguishable from eachother they averaged and kept. 
If the percentage difference is within 1% of the mean them, there are also kept. 
Otherwise their values are set to `NA`. 
} 
\arguments{ 
  \item{file_adiposity}{CrescNet file path for obese labeled subjects.} 
  \item{file_control}{CrescNet file path for control labeled subjects.} 
  \item{n_train}{Number of samples used for model building.} 
  \item{center_age}{Should age be centered? Default: FALSE} 
  \item{use_days}{Age given in days, instead of years? Default: FALSE} 
  \item{encode_sex}{Sex encoded as 1 (male) and -1 (female)? Default: FALSE} 
  \item{short_ids}{Shorter subject IDs? Default: TRUE} 
  \item{first_extra}{First measurements as separate variables? Default: FALSE} 
  \item{check_dups}{Check for duplicates? Default: FALSE} 
  \item{seed}{Seed used for random sampling of subjects. Default: 1} 
  \item{age_range}{Range of ages to keep. Default: NULL} 
} 
\value{ 
 A data.frame contraining all variables of the CrescNet dataset. 
} 
