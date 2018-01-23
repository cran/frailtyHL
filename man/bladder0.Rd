\name{bladder0}
\alias{bladder0}
\docType{data}
\title{Bladder cancer data}
\description{
Bladder0 is a subset of 410 patients from a full data set with bladder cancer
from 21 centers that participated in the EORTC trial
(Sylvester et al., 2006). Time to event is the duration of the disease free
interval (DFI), which is defined as time from randomization to the date of
the first recurrence.
}
\usage{data("bladder0")}
\format{
  A data frame with 410 observations on the following 5 variables.
  \describe{
    \item{\code{Center}}{Institution number of 24 centers}
    \item{\code{Surtime}}{Time to the first recurrence from randomization}
    \item{\code{Status}}{Censoring indicator(1=recurrence, 0=no event)}
    \item{\code{Chemo}}{Treatment indicator representing chemotherapy(0=No, 1=Yes)}
    \item{\code{Tustat}}{Indicator representing prior recurrent rate(0=Primary, 1=Recurrent)}
  }
}
\references{
Sylvester, R., van der Meijden, A.P.M., Oosterlinck, W., Witjes, J., Bouffioux, C., Denis, L., Newling, D.W.W. and Kurth, K. (2006).
Predicting recurrence and progression in individual patients with stage Ta T1 bladder cancer using EORTC risk tables: a combined analysis of
2596 patients from seven EORTC trials. European Urology, 49, 466-477.

Ha, I.D., Sylvester, R., Legrand, C. and MacKenzie, G. (2011).
Frailty modelling for survival data from multi-centre clinical trials. Statistics in Medicine, 30, 28-37.
}

