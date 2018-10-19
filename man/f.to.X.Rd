\name{f.to.X}
\alias{f.to.X}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data transformation: from species frequency counts to species abundance data
}
\description{
This function is to convert species frequency counts data to species abundance data.
}
\usage{
f.to.X(f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
Species frequency counts data.
}
}
\value{
Species abundance data is returned.
}
\author{
Youhua Chen & Tsung-Jen Shen
}

\seealso{
\code{\link{X.to.f}}
}

\examples{
## The sample is composed of 5 singletons, two doubletons, and one tripleton species.	
f = c(5, 2, 1)
f.to.X(f)
}
