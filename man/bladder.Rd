\name{bladder}
\alias{bladder}
\docType{data}
\title{Bladder Cancer Data}
\description{Bladder is an extension of Bladder0 to competing risks with 396 patients with
bladder cancer from 21 centers, focusing on two competing endpoints, i.e, time to
first bladder recurrence (an event of interest; Type 1 event) and time to death prior
to recurrence (competing event; Type 2 event).
}
\usage{data("bladder")}
\format{
  A data frame with 396 observations on the following 13 variables.
  \describe{
    \item{\code{OBS}}{Observation number}
    \item{\code{center}}{Institution number of 24 centers}
    \item{\code{surtime}}{Time to event}
    \item{\code{status}}{Event indicator(1=recurrence, 2=death before recurrence, 0=no event)}
    \item{\code{CHEMO}}{Treatment indicator representing chemotherapy(0=No, 1=Yes)}
    \item{\code{AGE}}{Age(0, <= 65 years; 1, > 65 years)}
    \item{\code{SEX}}{Sex(0=male, 1=female)}
    \item{\code{PRIORREC}}{Prior recurrent rate(0, primary; 1, <= 1/yr; 2, > 1/yr)}
    \item{\code{NOTUM}}{Number of tumors(0, single; 1, 2-7 tumors; 2, >= 8 tumors)}
    \item{\code{TUM3CM}}{Tumor size(0, < 3cm; 1, >= cm)}
    \item{\code{TLOCC}}{T cotegory(0=Ta, 1=T1)}
    \item{\code{CIS}}{Carcinoma in situ (0=No, 1=Yes)}
    \item{\code{GLOCAL}}{G grage(0=G1, 1=G2, 2=G3)}
  }
}
\references{
Sylvester, R., van der Meijden, A.P.M., Oosterlinck, W., Witjes, J., Bouffioux, C., Denis, L., Newling, D.W.W. and Kurth, K. (2006).
Predicting recurrence and progression in individual patients with stage Ta T1 bladder cancer using EORTC risk tables: a combined analysis of
2596 patients from seven EORTC trials. European Urology, 49, 466-477.

Ha, I.D., Sylvester, R., Legrand, C. and MacKenzie, G. (2011).
Frailty modelling for survival data from multi-centre clinical trials. Statistics in Medicine, 30, 28-37.
}
