\name{prepare_dataset} 
\alias{prepare_dataset} 
\title{Loading and preprocessing of the CrescNet dataset.} 
\usage{ 
prepare_dataset(file_adiposity, file_nonobese, n_train=400, is_znorm=TRUE, 
                center_age=FALSE, use_days=FALSE, encode_sex=FALSE, short_ids=TRUE, 
                first_extra=FALSE, check_dups=TRUE, seed=1234, set_age_range=TRUE, 
                body_status="IOTF", ...) 
} 
\arguments{ 
  \item{file_adiposity}{CrescNet file path for obese labeled subjects.} 
  \item{file_nonobese}{CrescNet file path for nonobese labeled subjects.} 
  \item{n_train}{Number of samples used for model building. Default: 400} 
  \item{is_znorm}{Include z-normalization? Default: TRUE} 
  \item{center_age}{Should age be centered? Default: FALSE} 
  \item{use_days}{Age given in days, instead of years? Default: FALSE} 
  \item{encode_sex}{Sex encoded as 1 (male) and -1 (female)? Default: FALSE} 
  \item{short_ids}{Shorter subject IDs? Default: TRUE} 
  \item{first_extra}{First measurements as separate variables? Default: FALSE} 
  \item{check_dups}{Check for duplicates? Default: TRUE} 
  \item{seed}{Seed used for random sampling of subjects. Default: 1234} 
  \item{set_age_range}{Should age be categorized in age categories? Default: TRUE} 
  \item{body_status}{The way of how BMI should be categorized. Possible values: "IOTF" or "adjusted". Default: "IOTF"} 
  \item{...}{Extra arguments to the function 'cutage'.} 
} 
\value{ 
 A list with the entries `train` and `test`, each containing a data.table with the CresNet data. 
} 
\description{Prepare the CrescNet dataset.} 
\details{ 
The function expects file paths to specific CrescNet data with the following variables: 
'subject_id', 'sex', 'gestational_age', 'age', 'height', 'weight', 'height_sds', 
'bmi', 'bmi_sds' and 'age_group'. `check_dups` will look for duplicated values. 
If they are indistinguishable from eachother, i.e. if the percentage difference is within 1% of their mean, 
there are averaged and kept. Otherwise they are set to be `NA`. 
 
The dataset will be splitted in a training and 
a test set by default, where `n_train` specifies the number of subjects in the training 
set, such that the return value is a list. By default, 400 subjects are assigned to the 
training set, while 1600 are assigned to the test test set. 
 
The argument `body_status` specifies, which kind of cutoffs should be used for defining 
body status based on BMI: The value `IOTF` uses the cutoffs based on the broad categories defined by 
the International Obesity Task Force (IOTF) guidelines: thinness (< 18.5 kg/m^2), 
normal (18.5-25.0 kg/m^2), overweight (>= 25.0 kg/m^2) and obesity (>= 30.0 kg/m^2), 
without finer distinction. Specifying the argument to `adjusted` applies the proposal of 
Onis and Lobstein (2010) to categorize body shape by BMI-SDS: 
Using +2 and +3 BMI-SDS for ages <= 5 years and +1 and +2 BMI-SDS thereafter for 
overweight and obesity, respectively. 
This is a consequence of the association of BMI to actual fat mass, which 
changes throughout the development of children. Thus, BMI cutoffs are not independent of age and 
using them indiscriminately, leads to overdiagnosing very small children as overweight / obese. 
Note, that this will result in a discontinuity, when the BMI body status for longitudinal 
data of individuals is tracked, which exemplifies the general problem of defining cutoffs for transforming 
continuous variables into ordinal categories (Harrell, 2015, p.18). 
} 
\references{ 
  \enumerate{ 
   \item{Onis, M. de and T. Lobstein (2010). “Deﬁning obesity risk status in the general childhood population: which cut-oﬀs should we use?” eng. In: International journal of pediatric obesity : IJPO : an oﬃcial journal of the International Association for the Study of Obesity. Vol. 5. England, pp. 458–60.} 
   \item{Harrell, Frank E. Jr. (Aug. 2015). Regression Modeling Strategies. With Applications to Linear Models, Logistic and Ordinal Regression, and Survival Analysis. 2nd ed. Springer Cham. isbn: 978-3-319-19424-0.} 
  } 
} 
\seealso{ 
\code{\link{cutage}}, \code{\link{znorm}}, \code{\link{linref}} 
} 
