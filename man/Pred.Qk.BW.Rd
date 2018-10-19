\name{Pred.Qk.BW}
\alias{Pred.Qk.BW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Incidence-based data: Bayesian-weight estimator
}
\description{
Bayesian-weight estimator for predicting the number of new rare species using incidence/quadrat data
}
\usage{
Pred.Qk.BW(Q, nT, u, b, k.show = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{
A vector of species frequency counts, i.e., the number of species dectected once (in only one quadrat), the number of species dectected twice (in exactly two quadrats), and so forth. 
}
  \item{nT}{
The number of quadrats of the original sample
}
  \item{u}{
The number of quadrats of an additional sample
}
  \item{b}{
A vector of two estimated parameters for obtaining the estimated relative species abundances by Chao et al.'s (2015) method.
}
  \item{k.show}{
Display the estimating results of the numbers of new rare species detected in the number of quadrats <= k.show in the additional sample
}
}
\value{
The numbers of new rare species detected in the number of quadrats <= k.show are estimated by the incidence-based Bayesian-weight estimator and returned. 
}
\references{
Shen TJ, Chen YH (2018) A Bayesian weighted approach to predicting the number of newly discovered rare species. Conservation Biology, In press.

}
\author{
Youhua Chen & Tsung-Jen Shen
}
\seealso{
\code{\link{Pred.Fk.BW}}
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
b = DetInc(y=Xi, nT=nT)	
Pred.Qk.BW(Q=Q, nT=nT, u=u, b=b[1:2])	
}

