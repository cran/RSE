\name{boot.incidence.fun}
\alias{boot.incidence.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a bootstrap incidence-based sample
}
\description{
Given an incidence-based data, a bootstrap sample is generated from a reconstructed bootstrap assemblage.
}
\usage{
boot.incidence.fun(S.hat, nT, Q, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S.hat}{
An estimate of species richness.
}
  \item{nT}{
The number of quadrats of the original sample
}
  \item{Q}{
A vector of species frequency counts, i.e., the number of species dectected once (in only one quadrat), the number of species dectected twice (in exactly two quadrats), and so forth. 
}
  \item{b}{
A vector of estimates of two parameters for obtaining the estimated detection probabilities of observed species in a given sample by Chao et al.'s (2015) method.
}
}
\value{
The generated bootstrap sample is a vector of species frequency counts, i.e., the number of species dectected once (in only one quadrat of the bootstrap sample), the number of species dectected twice (in exactly two quadrats of the bootstrap sample), and so forth. 
}
\references{
Chao A, Hsieh T, Chazdon R, Colwell R, Gotelli N. 2015. Unveiling the species-rank abundance distribution by generalizing the Good-Turing sample coverage theory. Ecology 96:1189-1201.
}
\author{
Youhua Chen & Tsung-Jen Shen
}
\seealso{
\code{\link{boot.abundance.fun}}
}
\examples{
## As an example, Canadian-mite data are used here.	
data(CanadaMite)
## two columns represent two samples of incidence counts
X.merge = CanadaMite
## the first column is treated as the sample
X.col1 = X.merge[,1]
Xi = X.col1
## Convert species incidence count data to frequency counts data
Q = X.to.f(Xi)
## the number of quadrats in the first sample
nT = 16
b = DetInc(Xi, nT)	
boot.incidence.fun(S.hat=sum(Q)+b[3], nT=nT, Q=Q, b=b[1:2])
}

