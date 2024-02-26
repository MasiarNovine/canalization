\name{do_simulation} 
\alias{do_simulation} 
\title{Run a Monte-Carlo sampling simulation to emulate CrescNet data selection} 
\description{ 
Investigation of the influence of sampling on the generated sample. 
} 
\details{ 
Simulation function to investigate the effect of the random stratified sampling, i.e. 
the use of split criteria based on centiles on the association of the BMI mean and standard deviation (SD). 
} 
\arguments{ 
\item{cresc_data:}{A `data.frame` with CrescNet data.} 
\item{n_subjects:}{Number of simulated subjects. Default: 1000.} 
\item{k:}{Number of observations per subject. Default: 20.} 
\item{n_pop:}{Number of observations in the population. Default: 1e5.} 
\item{split:}{Split criterion. Default: 0.97.} 
} 
\value{ 
A `list` with two data.tables, which include the simulated data `sample` and the summary statistics for each individual `aggregated`. 
} 
