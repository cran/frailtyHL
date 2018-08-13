\name{jointmodeling}
\alias{jointmodeling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Defining the Fixed and Random Models for the Mean and Dispersion Parameters in Joint Models}
\description{The jointmodeling specifies jointly both the hazard model in the frailty model and 
the mean model in the linear mixed model.}
\usage{
jointmodeling(Model = "mean", RespDist = "gaussian", Link = NULL, LinPred = "constant", 
RandDist = NULL, Offset = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Model}{This option specifies the mean model when Model="mean" (default).}
  \item{RespDist}{This option specifies the distribution of response variables (linear mixed model: "gaussian" 
  or accelerated failure time model : "AFT" or frailty model : "FM")}
  \item{Link}{The link function for the linear predictor is
specified by the option Link. For "AFT" or "FM" (or "gaussian") in RespDist, it is specified by "log" (or "identity").}
  \item{LinPred}{The option LinPred specifies the fixed and
random terms for the linear predictor.}
  \item{RandDist}{The option RandDist specifies the
distributions of the random terms represented in the option LinPred.}
  \item{Offset}{The option Offset can be used to
specify a known component to be included in the linear predictor specified by LinPred
during fitting.}
}
