\name{Pred.abundance.rare}
\alias{Pred.abundance.rare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Abundance-based data: predicting the number of new rare species
}
\description{
Abundance-based prediction on the number of new rare species using a Bayesian-weight and
two unweighted estimators along with their bootstrap standard errors and 95\% bootstrap confidence intervals.
}
\usage{
Pred.abundance.rare(boot.rep = 100, f = NULL, xi = NULL, m, k.show = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{boot.rep}{
Replicate number of the bootstrapping procedure
}
  \item{f}{
A vector of species frequency counts, i.e., the number of singleton species (only one individual observed in the sample), the number of doubleton species (two individuals observed in the sample), and so forth. 
}
  \item{xi}{
A vector of species abundance data, i.e., the number of individuals of species 1, the number of individuals of species 2, and so forth. 
}
  \item{m}{
The number of individuals of an additional sample
}
  \item{k.show}{
Display the estimating result of the numbers of extremely rare species with abundance <= k.show in the additional sample
}
}
\value{
Estimating results including point estimate, bootstrap standard error, and 95 \% bootstrap confidence interval for each of three methods (a Bayesian-weight and two unweighted estimators)
}
\references{
Shen TJ, Chen YH (2018) A Bayesian weighted approach to predicting the number of newly discovered rare species. Conservation Biology, In press.

}
\author{
Youhua Chen & Tsung-Jen Shen
}
\seealso{
\code{\link{Pred.incidence.rare}}
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
Pred.abundance.rare(f=f, m=m)	
}
