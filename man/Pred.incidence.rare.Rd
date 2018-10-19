\name{Pred.incidence.rare}
\alias{Pred.incidence.rare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Incidence-based data: predicting the number of new rare species
}
\description{
Incidence-based prediction on the number of new rare species using a Bayesian-weight and
two unweighted estimators along with their bootstrap standard errors and 95\% bootstrap confidence intervals.
}
\usage{
Pred.incidence.rare(boot.rep = 100, Q = NULL, xi = NULL, nT, u, k.show = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{boot.rep}{
Replicate number of the bootstrapping procedure
}
  \item{Q}{
A vector of species frequency counts, i.e., the number of species dectected once (in only one quadrat), the number of species dectected twice (in exactly two quadrats), and so forth. 
}
  \item{xi}{
A vector of species incidence counts, i.e., the number of quadrats with species 1, the number of quadrats with species 2, and so forth. 
}
  \item{nT}{
The number of quadrats of the original sample
}
  \item{u}{
The number of quadrats of an additional sample
}
  \item{k.show}{
Display the estimating results of the numbers of new rare species detected in the number of quadrats <= k.show in the additional sample
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
\code{\link{Pred.abundance.rare}}
}
\examples{
## As an example, Canadian-mite data are used here.	
data(CanadaMite)
## two columns represent two samples of incidence counts
X.merge = CanadaMite
## the first column is treated as the original sample
X.col1 = X.merge[,1]
Xi = X.col1
## Convert species incidence count data to frequency counts data
Q = X.to.f(Xi)
## the number of quadrats in the first sample
nT = 16
## the number of quadrats in the additional sample (i.e., the second column)
u = 16
Pred.incidence.rare(Q=Q, nT=nT, u=u)
}

