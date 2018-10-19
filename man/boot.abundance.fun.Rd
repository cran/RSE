\name{boot.abundance.fun}
\alias{boot.abundance.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a bootstrap abundance-based sample}
\description{
Given an abundance-based data, a bootstrap sample is generated from a reconstructed bootstrap assemblage.
}
\usage{
boot.abundance.fun(S.hat, f, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S.hat}{
An estimate of species richness.
}
  \item{f}{
A vector of species frequency counts, i.e., the number of singleton species (only one individual observed in the sample), the number of doubleton species (two individuals observed in the sample), and so forth. 
}
  \item{b}{
A vector of estimates of two parameters for obtaining the estimated relative abundances of observed species in a given sample by Chao et al.'s (2015) method.
}
}
\value{
The generated bootstrap sample is a vector of species frequency counts, i.e., the number of singleton species (only one individual observed in the bootstrap sample), the number of doubleton species (two individuals observed in the bootstrap sample), and so forth.
}
\references{
Chao A, Hsieh T, Chazdon R, Colwell R, Gotelli N. 2015. Unveiling the species-rank abundance distribution by generalizing the Good-Turing sample coverage theory. Ecology 96:1189-1201.
}
\author{
Youhua Chen & Tsung-Jen Shen
}
\seealso{
\code{\link{boot.incidence.fun}}
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
b = DetAbu(x=Xi, zero=FALSE)		
## the estimated number of species 
S.hat = SpEst.Chao1.abun(f)	
boot.abundance.fun(S.hat=S.hat, f=f, b=b)
}

