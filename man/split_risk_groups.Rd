\name{split_risk_groups} 
\alias{split_risk_groups} 
\title{Split the at-risk group based on percentiles} 
\description{ 
Split the based on given cutoff values. 
} 
\details{ 
Extract respective subgroups based on given percentiles for 
BMI-SDS values. For each subject, it is checked, whether its endpoint 
measurement is >97th percentile or not. Additionally, it is checked, 
whether the last measurement is <=85th percentile, i.e. whether individual 
subjects are non-obese concerning their BMI-SDS. The result is a list with 
two `data.frames`: 
One with the extracted values for the at-risk (non-obese) and one 
for the at-risk (obese) group. The idea is to have a projective outlook 
on the samples based on the respective 
end measurement. Note, that the end measurement for each subject might 
not be the same, but it lies approximately in the range of 14 to 18 years. 
The values can be plotted to observe the change of BMI-SDS in the 
course of the development of each subject. 
} 
\arguments{ 
  \item{cresc_data}{A `data.frame` with CrescNet data.} 
  \item{use_p85}{Should the upper bound for the non-obese group be the 85th percentile? Default: TRUE} 
} 
\value{ 
  A `list` with two objects of class `data.frame`. 
} 
