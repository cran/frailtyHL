\name{cgd}
\Rdversion{1.0}
\alias{cgd}
%\alias{frailtyHL}

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
    \item{\code{hospital}}{Hospital number for 13 hospitals.}
    \item{\code{patient}}{Patient number for 128 patients.}
    \item{\code{time}}{The recurrent infection times of each pateient from the different hospitals.}
    \item{\code{censor}}{Censoring indicator (1=uncensored, 0=censored).}
    \item{\code{trt}}{Treatment (1=gamma interferon, 0=placebo).}
    \item{\code{inherit}}{Pattern of inheritance (0=autosomal recessive, 1=X-linked).}
    \item{\code{age}}{Age of each patient, in years.}
    \item{\code{height}}{Height of each patient, in cm.}
    \item{\code{weight}}{Weight of each patient, in kg.}
    \item{\code{steroids}}{Using corticosteroids at times of study centry(1=Yes, 0=No).}
    \item{\code{proylac}}{Using prophylactic antibiotics at time of study entry(1=Yes, 0=No).}    
    \item{\code{sex}}{Sex of each patient(1=male, 0=female).}    
    \item{\code{h}}{Hospital region(0=U.S.A., 1=Europe).}    
    \item{\code{longi}}{A longitudinal variable representing the accumulated time from the first infection(in years).}        
    \item{\code{enum}}{Sequence number. For each patient, the infection records are in sequnce number order.}            
  }
}


\references{Fleming, T. R. and Harrington, D. R. (1991). Counting processes and survival analysis. Wiley: New York.
}
\examples{
data(cgd)
}

