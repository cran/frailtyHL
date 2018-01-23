\name{renal}
\alias{renal}
\docType{data}
\title{Renal transplant data}
\description{
This is a data set from a clinical study to investigate the chronic renal allograft
dysfunction in renal transplants (Ha et al., 2018).
Data were available from 87 male and 25 female renal transplanted patients
who survived more than 4 years after transplant. For each patient,
both repeated-measure outcomes (serum creatinine levels) at several time points
and a terminating event time (graft-loss time) were observed.
}
\usage{data("renal")}
\format{
  A data frame with 1395 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{Patient id}
    \item{\code{month}}{Time points (month) at which the measurements of sCr were recorded}
    \item{\code{cr}}{Serum creatinine (sCr) level}
    \item{\code{sex}}{Sex(1=male, 0=female)}
    \item{\code{age}}{Age(years)}
    \item{\code{icr}}{Reciprocal of sCr(=1/sCr)}
    \item{\code{sur_time}}{Time to graft loss}
    \item{\code{status}}{Censoring indicator(1=graft loss, 0=no event)}
    \item{\code{first}}{The first survival time (time to graft loss) of each patient}
  }
}
\references{
Ha, I. D., Noh, M. and Lee, Y. (2017). H-likelihood approach for joint modelling of longitudinal outcomes and time-to-event data.
Biometrical Journal, 59, 1122--1143.

Ha, I. D., Jeong, J.H. and Lee, Y. (2018). Statistical modelling of survival data with
random effects: h-likelihood approach. Springer.
}

