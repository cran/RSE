\name{Pred.Fk.unweighted}
\alias{Pred.Fk.unweighted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Abundance-based data: Unweighted estimator
}
\description{
Unweighted estimator based on Chao et al. (2015)'s paper using abundance-based data for predicting the number of new rare species in an additional ecological sample
}
\usage{
Pred.Fk.unweighted(f, m, b, f0, k.show = 3)
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
  \item{f0}{
The estimated number of unseen species in the original sample by Chao 1 estimator (Chao 1984)
}
  \item{k.show}{
Display the estimating result of the numbers of extremely rare species with abundance <= k.show in the additional sample
}
}
\value{
The numbers of new rare species with abundance <= k.show are estimated by the abundance-based unweighted estimator and returned. 
}
\references{
Chao A, Hsieh T, Chazdon R, Colwell R, Gotelli N. 2015. Unveiling the species-rank abundance distribution by generalizing the Good-Turing sample coverage theory. Ecology 96:1189-1201.
	
Chao A. 1984. Non-parametric estimation of the number of classes in a population. Scandinavian Journal of Statistics 11:265-270.

Shen TJ, Chen YH (2018) A Bayesian weighted approach to predicting the number of newly discovered rare species. Conservation Biology, In press.

}
\author{
Youhua Chen & Tsung-Jen Shen
}
\seealso{
\code{\link{Pred.Qk.unweighted}}
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
## the estimated number of unseen species in the original sample
f0 = SpEst.Chao1.abun(f)-sum(f)	
Pred.Fk.unweighted(f=f, m=m, b=b, f0=f0)	
}

