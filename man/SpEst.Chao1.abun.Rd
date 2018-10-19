\name{SpEst.Chao1.abun}
\alias{SpEst.Chao1.abun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Species richness estimation
}
\description{
Chao1 estimator of species richness
}
\usage{
SpEst.Chao1.abun(f)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
A vector of species frequency counts, i.e., the number of singleton species (only one individual observed in the sample), the number of doubleton species (two individuals observed in the sample), and so forth. 
}
}
\references{
Chao A. 1984. Non-parameteric estimation of the number of classes in a population. Scandinavian Journal of Statistics 11:265-270
}
\author{
Youhua Chen & Tsung-Jen Shen
}

\examples{
## As an example, Herpetological assemblage data are used here.		
data(HerpetologicalData)
## two columns represent two samples of species abundance data
X.merge = HerpetologicalData
## the first column is treated as the original sample
Xi = X.merge[,1]
## Convert species abundance data to species frequency counts data
f = X.to.f(Xi)
SpEst.Chao1.abun(f=f)
}

