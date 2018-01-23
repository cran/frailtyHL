\name{frailtyHL-package}
\alias{frailtyHL-package}
\docType{package}
\title{H-likelihood Approach for Frailty Models}
\description{The frailtyHL package fits frailty models which are Cox's proportional hazards models incorporating
random effects. The function implements the h-likelihood estimation procedures. For the frailty
distribution lognormal and gamma are allowed. The h-likelihood uses the Laplace approximation
when the numerical integration is intractable, giving a statistically efficient estimation in frailty
models. (Ha, Lee and Song, 2001; Ha and Lee, 2003, 2005; Lee, Nelder and Pawitan, 2017; Ha, Jeong and Lee, 2018).
This package handles various random-effect survival models such as multi-lever or time-dependent frailties, 
competing-risk frailty models, accelerated failure time (AFT) random-effect models, and joint modelling of linear mixed models and frailty models.
It also provides penalized variable-selection procedures (LASSO, SCAD and HL).}
\details{
\tabular{ll}{
Package: \tab frailtyHL\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2018-01-02\cr
License: \tab Unlimited\cr
LazyLoad: \tab yes\cr
}
This is version 2.1 of the frailtyHL package.
}
\author{Il Do Ha, Maengseok Noh, Jiwoong Kim, Youngjo Lee

Maintainer: Maengseok Noh <msnoh@pknu.ac.kr>
}
\references{
Ha, I.D. and Lee, Y. (2003). Estimating frailty models via Poisson Hierarchical generalized linear
models. Journal of Computational and Graphical Statistics, 12, 663-681.

Ha, I.D. and Lee, Y. (2005). Comparison of hierarchical likelihood versus orthodox best linear
unbiased predictor approaches for frailty models. Biometrika, 92, 717-723.

Ha, I.D., Lee, Y. and Song, J. K. (2001). Hierarchical likelihood approach for frailty models.
Biometrika, 88, 233-243.

Ha, I.D., Jeong, J.H. and Lee, Y. (2018). Statistical modelling of survival data with random effects.
Springer.

Lee, Y., Nelder, J.A. and Pawitan, Y. (2017). Generalised linear models with random effects:
unified analysis via h-likelihood. 2nd Edition. Chapman & Hall: London.}

