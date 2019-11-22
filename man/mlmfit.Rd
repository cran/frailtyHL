\name{mlmfit}
\alias{mlmfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accelerated Failure Time (AFT) Models with Random Effects}
\description{mlmfit is used to fit linear mixed models with censoring 
by using h-likelihood.}
\usage{
mlmfit(jm1, data, weights, subset, na.action, Maxiter = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{jm1}{This option requires \code{jointmodeling} object which specifies the AFT random-effect model.}
  \item{data}{dataframe containing the variables used in the jm1}
    \item{weights}{Vector of case weights.}
  \item{subset}{Expression indicating which subset of the rows of data should be used in the fit.
All observations are included by default.}
  \item{na.action}{A missing-data filter function.}
  \item{Maxiter}{numeric, maximum number of iterations}
}

