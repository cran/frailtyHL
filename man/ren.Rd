\name{ren}
\alias{ren}
\docType{data}
\title{Mammary tumor data}
\description{
The data set by presented Gail et al. (1980) is based on multiple occurrences of mammary tumors for 48
female rats. The primary outcome of interest was time to development of a mammary
tumor for 23 female rats in the treatment group and 25 female rats in the control
group. Initially, 76 rats were injected with a carcinogen for mammary cancer at
day zero, and then all rats were given retinyl acetate to prevent cancer for 60 days.
After 60 days, forty-eight rats which remained tumor-free were randomly assigned to
continue being treated with retinoid prophylaxis (treatment group) or to the control
group receiving no further retinoid prophylaxis. Rats were palpated for tumors twice
weekly and observation ended 182 days after the initial carcinogen injection.
In some cases, there were multiple tumors detected by the same day.
The number of tumors ranges from 0 to 13.

}
\usage{data("ren")}
\format{
  A data frame with 254 observations on the following 6 variables.
  \describe{
    \item{\code{rat}}{Rat id}
    \item{\code{time1}}{Start time}
    \item{\code{time2}}{Stop time}
    \item{\code{del}}{Censoring indicator(1=tumor, 0=censored)}
    \item{\code{gp}}{Treatment indicator(1=drug, 0=control)}
    \item{\code{time}}{time2-time1 (time=time+0.01 if there are ties)}
  }
}
\references{
Gail, M.H. Santner, T.J. and Brown, C.C. (1980), An analysis of comparative carcinogenesis experiments
based on multiple times to tumor. Biometrics, 36, 255-266.

Ha, I. D., Jeong, J.H. and Lee, Y. (2018). Statistical modelling of survival data with
random effects: h-likelihood approach. Springer.
}
