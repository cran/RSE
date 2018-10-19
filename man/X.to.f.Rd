\name{X.to.f}
\alias{X.to.f}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data transformation: from species abundance data to species frequency counts data
}
\description{
This function is to convert a vector of species abundance data to a vector of species frequency counts data.
}
\usage{
X.to.f(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A vector of species abundance data.
}
}

\value{
Species frequency counts is returned.
}
\author{
Youhua Chen & Tsung-Jen Shen
}
\seealso{
\code{\link{f.to.X}}
}
\examples{
## As an example, Herpetological assemblage data are used here.		
data(HerpetologicalData)
## two columns represent two samples of species abundance data
X.merge = HerpetologicalData
## the first column is treated as the original sample
X.col1 = X.merge[,1]
Xi = X.col1
## convert species abundance data to species frequency counts data
X.to.f(Xi)	
}

