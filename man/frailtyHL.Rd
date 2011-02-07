\name{frailtyHL}
\Rdversion{1.1}
\alias{frailtyHL}

\docType{package}

\title{Fitting Frailty Models using H-likelihood Approach}

\description{frailtyHL is used to fit frailty models using h-likelihood estimation procedures. 
For the frailty distribution lognormal and gamma are allowed. 
In particular, nested (multilevel) frailty models allow survival studies for hierarchically clustered data by including two iid normal random effects.
The h-likelihood uses the Laplace approximation when the numerical integration is intractable, giving a statistically efficient estimation in frailty models.
(Ha, Lee and Song, 2001; Ha and Lee, 2003, 2005; Lee, Nelder and Pawitan, 2006).}

\usage{
frailtyHL(formulaMain,censor,DataMain,RandDist="Normal",mord=0,dord=1,
Maxiter=200,convergence=10^-6)
}

\arguments{
  \item{formulaMain}{Specify the formula for the mean structure of the model. e.g. y ~ x + (1|id), y : response, x : fixed covariate, id : random effect.}
  \item{censor}{One vector for censoring indicator (1=uncensored, 0=censored).}
  \item{DataMain}{Dataframe for formulaMain.}
  \item{RandDist}{Distribution for random effect ("Normal" or "Gamma").}
  \item{mord}{The order of Laplce approximation to fit the mean parameters (0 or 1); default=0.}
  \item{dord}{The order of Laplace approximation to fit the dispersion components (1 or 2); default=1.}
  \item{Maxiter}{The maximum number of iterations; default=200.}
  \item{convergence}{Specify the convergence criterion, the default is 1e-6.}
}

\details{
frailtyHL package produces estimates of fixed effects and frailty parameters as well as their standard errors. 
Also, frailtyHL makes it possible to fit models where the frailty distribution is normal and gamma and estimate variance components 
when frailty structure is allowed to be shared or nested.
}

\references{
Ha, I. D. and Lee, Y. (2003). Estimating frailty models via Poisson Hierarchical generalized linear models.
Journal of Computational and Graphical Statistics, 12, 663--681.

Ha, I. D. and Lee, Y. (2005). Comparison of hierarchical likelihood versus orthodox best linear unbiased predictor approaches
for frailty models. Biometrika, 92, 717--723.

Ha, I. D., Lee, Y. and Song, J. K. (2001). Hierarchical likelihood approach for frailty models. Biometrika, 88, 233--243.

Lee, Y., Nelder, J. A. and Pawitan, Y. (2006). Generalised linear models with random effects: unified analysis via
h-likelihood. Chapman & Hall: London.
}
\keyword{ package }
\seealso{
<\code{\link{summary.frailtyHL}}>
}

\examples{
#### Analysis of Kidney data
data(kidney)
#### Normal frailty model using order = 1, 1 for the mean and dispersion
kidney_ln11<-frailtyHL(time~sex+age+(1|patient)-1,kidney$censor,kidney,
RandDist="Normal",mord=1,dord=1)
#### Normal frailty model using order = 1, 2 for the mean and dispersion
kidney_ln12<-frailtyHL(time~sex+age+(1|patient)-1,kidney$censor,kidney,
RandDist="Normal",mord=1,dord=2)
#### Gamma frailty model using order = 1, 1 for the mean and dispersion
kidney_g11<-frailtyHL(time~sex+age+(1|patient)-1,kidney$censor,kidney,
RandDist="Gamma",mord=1,dord=1)
#### Gamma frailty model using order = 1, 2 for the mean and dispersion
kidney_g12<-frailtyHL(time~sex+age+(1|patient)-1,kidney$censor,kidney,
RandDist="Gamma",mord=1,dord=2)

#### Analysis of CGD data
data(cgd)
#### Multilevel normal frailty model using order = 1, 1 for the mean and dispersion
cgd11<-frailtyHL(time~trt+(1|hospital)+(1|patient)-1,cgd$censor,cgd,
RandDist="Normal",mord=1,dord=1)
#### Multilevel normal frailty model using order = 1, 2 for the mean and dispersion
cgd12<-frailtyHL(time~trt+(1|hospital)+(1|patient)-1,cgd$censor,cgd,
RandDist="Normal",mord=1,dord=2)
}

