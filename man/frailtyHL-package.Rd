\name{frailtyHL-package}
\Rdversion{1.0}
\alias{frailtyHL-package}
\docType{package}
\title{H-likelihood Approach for Frailty Models}
\description{The frailtyHL package fits frailty models which are Cox's proportional hazards models incorporating random effects. 
The function implements the h-likelihood estimation procedures.
For the frailty distribution lognormal and gamma are allowed. 
The h-likelihood uses the Laplace approximation when the numerical integration is intractable, giving a statistically efficient estimation in frailty models.
(Ha, Lee and Song, 2001; Ha and Lee, 2003, 2005; Lee, Nelder and Pawitan, 2006).}
\details{
\tabular{ll}{
Package: \tab frailtyHL\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-02-07\cr
License: \tab Unlimited\cr
LazyLoad: \tab yes\cr
}
This is version 1.0 of the frailtyHL package.
}
\author{
Il Do Ha, Manegseok Noh, Youngjo Lee

Maintainer: Maengseok Noh <msnoh@pknu.ac.kr>
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
<\code{\link{frailtyHL}}>
}
\examples{
data(kidney)
kidney_g12<-frailtyHL(time~sex+age+(1|patient)-1,kidney$censor,kidney,
RandDist="Gamma",mord=1,dord=2)
}

