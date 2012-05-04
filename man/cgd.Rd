\name{cgd}
\Rdversion{1.1}
\alias{cgd}

\docType{package}

\title{Chronic Granulomatous Disease (CGD) Infection Data}

\description{The CGD data set in Fleming and Harrington (1991) is from a placebo-controlled randomized trial of gamma interferon 
in chronic granulomatous disease. In total, 128 patients from 13 hospitals were followed for about 1 year. 
The number of patients per hospital ranged from 4 to 26. Each patient may experience more than one infection. 
The survival times (times-to-event) are the times between recurrent CGD infections on each patient (i.e. gap times).  
Censoring occurred at the last observation for all patients, except one, who experienced a serious infection on the date he left the study.
}
\usage{data(cgd)}

\format{
  CGD data set contains 15 columns and 203 rows. A brief description of the data column is given below.
  \describe{
    \item{\code{id}}{Patient number for 128 patients.}
    \item{\code{center}}{Enrolling center number for 13 hospitals.}
    \item{\code{random}}{Date of randomization.}
    \item{\code{treat}}{Gamma-interferon treatment(rIFN-g) or placebo(Placebo).}
    \item{\code{sex}}{Sex of each patient(male, female).}    
    \item{\code{age}}{Age of each patient at study entry, in years.}
    \item{\code{height}}{Height of each patient at study entry, in cm.}
    \item{\code{weight}}{Weight of each patient at study entry, in kg.}
    \item{\code{inherit}}{Pattern of inheritance (autosomal recessive, X-linked).}
    \item{\code{steroids}}{Using corticosteroids at times of study centry(1=Yes, 0=No).}
    \item{\code{proylac}}{Using prophylactic antibiotics at time of study entry(1=Yes, 0=No).}    
    \item{\code{hos.cat}}{A categorization of the hospital region into 4 groups.}    
    \item{\code{tstart}}{Start of each time interval.}
    \item{\code{enum}}{Sequence number. For each patient, the infection records are in sequnce number order.}            
    \item{\code{tstop}}{End of each time interval.}
    \item{\code{status}}{Censoring indicator (1=uncensored, 0=censored).}
  }
}
\references{
Fleming, T. R. and Harrington, D. R. (1991). Counting processes and survival analysis. Wiley: New York.

Therneasu, T. (2012). survival: survival analysis, including penalised likelihood. http://CRAN.R-project.org/package=survival. R pakcage version 2.36-14.
}
\examples{
data(cgd)
}

