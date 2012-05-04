\name{summary.frailtyHL}
\Rdversion{1.1}
\alias{summary.frailtyHL}

\docType{package}

\title{Summary Method for frailtyHL Object}

\description{It provides standard summary statistics for the fitted frailtyHL objects.}

\usage{
\method{summary}{frailtyHL}(object, \dots)
}

\arguments{
  \item{object}{A \code{\link{frailtyHL}} object.}
  \item{...}{other arguments}
}

\details{
When applied to a frailtyHL object summary statistics are producted if the underlying model estimation converged. Otherwise, produces an error message.
}

\value{
  It returns a list of summary.frailtyHL object containing the following statistics and tables.
  \item{Model}{A selected model of log-normal or gamma frailty models.}
  \item{formula}{The formula used in the undelying \code{\link{frailtyHL}}.}
  \item{Method}{The order of Laplace approximation mord and dord used in the undelying \code{\link{frailtyHL}}.}
  \item{FixCoef}{A two-dimensional matrix of fixed effects parameter estimates, repective standard errors and t-values.}
  \item{RandCoef}{A two-dimensional matrix of random effects estimates and repective standard errors.}
  \item{likelihood}{Vectors of -2 * likleihood.}
  \item{iter}{Number of iterations used in \code{\link{frailtyHL}} estimation.}
  \item{convergence}{A text, "converged" if the estimation converged and "did not converge" oterwise.}
  \item{aic}{Three AIC criteria(cAIC, mAIC, rAIC) based on conditional likelihood, marginal likelihood and restricted likelihood}
}

