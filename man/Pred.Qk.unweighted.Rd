\name{Pred.Qk.unweighted}
\alias{Pred.Qk.unweighted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Incidence-based data: Unweighted Estimator}
\description{
Unweighted Estimator derived from Chao et al. (2015)'s paper using incidence/quadrat data for predicting the number of new rare species in an additional ecological sample
}
\usage{
Pred.Qk.unweighted(Q, nT, u, b, Q0, k.show = 3)
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
  \item{Q0}{
The estimated number of unseen species in the original sample by Chao 2 estimator (Chao 1987)
}
  \item{k.show}{
Display the estimating results of the numbers of new rare species detected in the number of quadrats <= k.show in the additional sample
}
}
\value{
The numbers of new rare species detected in the number of quadrats <= k.show are estimated by the incidence-based unweighted estimator derived from Chao et al. (2015)'s paper and returned. 
}
\references{
Chao A, Hsieh T, Chazdon R, Colwell R, Gotelli N. 2015. Unveiling the species-rank abundance distribution by generalizing the Good-Turing sample coverage theory. Ecology 96:1189-1201.

Chao A. 1987. Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43:783-791.

Shen TJ, Chen YH (2018) A Bayesian weighted approach to predicting the number of newly discovered rare species. Conservation Biology, In press.
}
\author{
Youhua Chen & Tsung-Jen Shen
}
\seealso{
\code{\link{Pred.Fk.unweighted}}
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
b = DetInc(Xi, nT)			
Pred.Qk.unweighted(Q=Q, nT=nT, u=u, b=b[1:2], Q0=b[3])		
}

