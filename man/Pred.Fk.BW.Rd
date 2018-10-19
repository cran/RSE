\name{Pred.Fk.BW}
\alias{Pred.Fk.BW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Abundance-based data: Bayesian-weight estimator
}
\description{
Bayesian-weight estimator for predicting the number of new rare species using abundance data as inputs
}
\usage{
Pred.Fk.BW(f, m, b, k.show = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{
A vector of species frequency counts, i.e., the number of singleton species (only one individual observed in the sample), the number of doubleton species (two individuals observed in the sample), and so forth. 
}
  \item{m}{
The number of individuals of an additional sample
}
  \item{b}{
A vector of two estimated parameters for obtaining the estimated relative species abundances by Chao et al.'s (2015) method.
}
  \item{k.show}{
Display the estimating result of the numbers of extremely rare species with abundance <= k.show in the additional sample
}
}
\value{
The numbers of new rare species with abundance <= k.show are estimated by the abundance-based Bayesian-weight estimator and returned. 
}
\references{
Shen TJ, Chen YH (2018) A Bayesian weighted approach to predicting the number of newly discovered rare species. Conservation Biology, In press.

}
\author{
Youhua Chen & Tsung-Jen Shen
}
\seealso{
\code{\link{Pred.Qk.BW}}
}
\examples{
## As an example, Herpetological assemblage data are used here.		
data(HerpetologicalData)
## two columns represent two samples of species abundance data
X.merge = HerpetologicalData
## the first column is treated as the original sample
X.col1 = X.merge[,1]
## the second column is treated as the additional sample
X.col2 = X.merge[,2]
Xi = X.col1
## Convert species abundance data to species frequency counts data
f = X.to.f(Xi)
## the number of individuals of the additional sample 
m = sum(X.col2)
b = DetAbu(x=Xi, zero=FALSE)		
Pred.Fk.BW(f=f, m=m, b=b)	
}

