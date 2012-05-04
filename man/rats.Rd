\name{rats}
\Rdversion{1.1}
\alias{rats}

\docType{package}

\title{Rats data}

\description{Rats data set presented by Mantel et al. (1977) is based on a tumorigenesis study of 50 litters of female rats.
For each litter, one rat was selected to receive the drug and the other two rats were placebo-treated controls. The survival
time is the time to the development of tumor, measured in weeks. Death before occurrent of tumor yields a right-censored obervation;
40 rats developed a tumor, leading to censoring of about 73 percent.
}
\usage{data(rats)}

\format{
  Rats data set contains 4 columns and 150 rows. A brief description of the data column is given below.
  \describe{
    \item{\code{litter}}{Litter number for 50 female rats.}
    \item{\code{rx}}{Treatment(1=drug, 0=placebo)}
    \item{\code{time}}{Time to the devlopment of tumor in weeks.}
    \item{\code{status}}{Censoring indicator(1=uncensored, 0=censored).}
  }
}


\references{
Mantel,N., Bohidar N. R. and Ciminera, J. L. (1977). Mantel-Haenszel analyses of litter-matched time-to-response data, 
with modifications for recovery of interlitter information. Cancer Research, 37, 3863--3868.

Therneasu, T. (2012). survival: survival analysis, including penalised likelihood. http://CRAN.R-project.org/package=survival. R pakcage version 2.36-14.
}
\examples{
data(rats)
}
