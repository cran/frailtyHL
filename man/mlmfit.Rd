\name{mlmfit}
\alias{mlmfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accelerated Failure Time (AFT) Models with Random Effects}
\description{mlmfit is used to fit linear mixed models with censoring 
(i.e. AFT models) by using h-likelihood.}
\usage{
mlmfit(jm1, data, Maxiter = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{jm1}{This option requires \code{jointmodeling} object which specifies the AFT random-effect model.}
  \item{data}{dataframe containing the variables used in the jm1}
  \item{Maxiter}{numeric, maximum number of iterations}
}

