\name{check_mixed_model} 
\alias{check_mixed_model} 
\title{Checking mixed-effects model} 
\description{Model check applying common checks for model assessment.} 
\usage{check_mixed_model(model, cresc_data, seed)} 
\arguments{ 
  \item{model}{A `gls` model.} 
  \item{cresc_data}{CrescNet data.table} 
  \item{seed}{The seed for extracting the 20 randomly drawn subjects. Default: 1} 
} 
\details{ 
The usually model checks are applied, i.e. (1) constant variance residuals check, (2) random distributed residuals, 
(3) goodness of fit by plotting fitted values against the outcome variable, (4) the variogram for observing the intra-subject correlation 
and (5) plotting fitted values against the standardized residuals. 
} 
\value{Series of `ggplot` plots.} 
