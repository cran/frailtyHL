\name{kidney}
\Rdversion{1.1}
\alias{kidney}

\docType{package}

\title{Kidney Infection Data}

\description{The data presented by McGilchrist and Aisbett (1991) consist of times to the first and second recurrences of infection 
in 38 kidney patients using a portable dialysis machine. Infections can occur at the location of insertion of the catheter. 
The catheter is later removed if infection occurs and can be removed for other reasons, in which case the observation is censored.
}
\usage{data(kidney)}

\format{
  Kidney data set contains 8 columns and 76 rows. A brief description of the data column is given below.
  \describe{
    \item{\code{id}}{Patient number for 38 patients.}
    \item{\code{time}}{Time to infection since insertion of the catheter}
    \item{\code{status}}{Censoring indicator(1=uncensored, 0=censored).}
    \item{\code{age}}{Age of each patient, in years.}
    \item{\code{sex}}{Sex of each patient(1=male, 2=female).}
    \item{\code{disease}}{Disease type(GN, AN, PKD, Other).}
    \item{\code{frail}}{Frailty estimate from original paper.}
  }
}


\references{
McGilchrist, C. A. and Aisbett, C. W. (1991). Regression with frailty in survival analysis. Biometrics, 47, 461--466.

Therneasu, T. (2012). survival: survival analysis, including penalised likelihood. http://CRAN.R-project.org/package=survival. R pakcage version 2.36-14.
}
\examples{
data(kidney)
}

