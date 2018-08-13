\name{CmpRsk}
\alias{CmpRsk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model Formula of Competing Risk }
\description{
A CmpRsk object is used as the response variable in the model formula. It is created using
the function CmpRsk(time, index), where time is the event time and index is an event indicator.
}
\usage{
CmpRsk(time, index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{the event time}
  \item{index}{the event indicator;
values of index must be sequential whole numbers where 0 denotes right censoring
and positive numbers refer to different event types.}
}

