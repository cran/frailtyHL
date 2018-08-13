\name{frailty.vs}
\alias{frailty.vs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Penalized Variable Selection for Frailty Models}
\description{frailty.vs is variable-selection procedures 
(LASSO, SCAD and HL) of fixed effects in frailty models.}
\usage{
frailty.vs(formula, model, penalty, data, B = NULL, v = NULL, 
alpha = NULL, tun1 = NULL, tun2 = NULL, varfixed = FALSE, varinit = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula object, with the response on the left of a ~ operator, and the terms for
the fixed and random effects on the right. e.g. formula=Surv(time,status)~x+(1|id),
time : survival time, status : censoring indicator having 1 (0) for uncensored
(censored) observation, x : fixed covariate, id : random effect.}
  \item{model}{Log-normal frailty models ("lognorm")}
  \item{penalty}{Penalty functions ("LASSO" or "SCAD" or "HL"))}
  \item{data}{Dataframe used}
  \item{B}{Initial values of fixed effects}
  \item{v}{Initial values of random effects. Zeros are default}
  \item{alpha}{Initial value of variance of random effects.}
  \item{tun1}{Tuning parameter gamma for LASSO, SCAD and HL}
  \item{tun2}{Tuning parameter omega for HL}
  \item{varfixed}{Logical value: if TRUE (FALSE), the value of one or more of the variance terms for the frailties is fixed (estimated).}
  \item{varinit}{Starting values for frailties, the default is 0.1.} 
}

