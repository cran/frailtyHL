\name{test}
\alias{test}
\docType{data}
\title{Simulated data with clustered competing risks}
\description{
A data set for the cause-specific hazard frailty model assuming a bivariate normal
distribution is generated using a technique similar to Beyersmann et al. (2009) and
Christian et al. (2016). Let there be two event types, Types 1 and 2, as well as
independent censoring. Consider a sample size n = 100 with (q, ni) = (50, 3).
Here, q is the number of clusters and ni is the cluster size.
The random effects (log-frailties) are from bivariate normal with mean vector (0,0) and variance-covariance matrix having (1,1,-0.5).
Data are generated from the conditional cause-specific hazard rates for each event type given the random effects.
Here, for Type 1 event the two true regression parameters are (0.6, -0.4) with a constant baseline hazard 2 and for Type 2 event
the true parameters are (-0.3, 0.7) with a constant baseline hazard 0.5, respectively.
The covariates x1 and x2 are generated from a standard normal distribution and a Bernoulli distribution  with probability 0.5, respectively.
Censoring times are generated from a Uniform(0, 1.3) distribution. Under this scenario,
with 25.2\% censoring, the proportions of Type 1 and Type 2 events are 53.2\% and 21.6\%, respectively.
}
\usage{data("test")}
\format{
  A data frame with 250 observations on the following 6 variables.
  \describe{
    \item{\code{obs}}{Observation number}
    \item{\code{id}}{Id number}
    \item{\code{time}}{Time to event}
    \item{\code{status}}{Event indicator(2=Type 2 event, 1=Type 1 event, 0=censored)}
    \item{\code{x1}}{A covariate from standard normal distribution}
    \item{\code{x2}}{A covariate from Bernoulli normal distribution}
  }
}
\references{
Beyersmann, J., Dettenkofer, M., Bertz, H. and Schumacher, M. (2007). A competing risks analysis of bloodstream infection after stem-cell
transplantation using subdistribution hazardsa and cause-specific hazards. Statistics in Medicine, 26, 5360-5369.

Christian, N. J., Ha, I. D. and Jeong, J. H. (2016). Hierarchical likelihood inference on clustered competing risks data. Statistics in Medicine, 35, 251-267.

Ha, I. D., Jeong, J. H. and Lee, Y. (2017). Statistical modelling of survival data with
random effects: h-likelihood approach. Springer, in press.
}

