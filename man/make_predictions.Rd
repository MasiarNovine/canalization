\name{make_predictions} 
\alias{make_predictions} 
\title{Generate predictions} 
\description{Uses the given GLS model to make predictions for 20 subjects.} 
\usage{make_predictions(model, cresc_data, seed, log)} 
\arguments{ 
  \item{model}{A `gls` model.} 
  \item{cresc_data}{CrescNet data.table} 
  \item{seed}{The seed for extracting the 20 randomly drawn subjects. Default: 1} 
  \item{log}{Whether the logarithm for the data should be used. Default: TRUE} 
} 
\details{ 
The function can be used to observe, how the model fits body weights of 20 randomely selected subjects of a test set. 
Next to a `ggplot` plot with fits for each subject, the estimated weight values are appended to the initilal data.table object. 
} 
\value{A `ggplot` object and a `data.table` with prediction estimates.} 
