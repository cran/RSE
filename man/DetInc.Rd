\name{DetInc}
\alias{DetInc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Incidence-based data: the estimation of parameters for obtaining the estimation of detection probabilites of observed species 
}
\description{
The estimation of parameters for obtaining the estimation of detection probabilites of observed species 
}
\usage{
DetInc(y, nT, zero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A vector of species incidence counts, i.e., the number of quadrats with species 1, the number of quadrats with species 2, and so forth. 
}
  \item{nT}{
The number of quadrats of the original sample
}
  \item{zero}{
A logical value for whether reserving zero frequency or not.
}
}
\value{
A vector of 3 elements: the first two values are the estimates of two parameters in Chao et al. (2015) for jointly estimating detection probabilities of observed species and the third one is the estimated number of unseen species in the sample by Chao 2 estimator (Chao, 1987).
}
\references{
Chao A, Hsieh T, Chazdon R, Colwell R, Gotelli N. 2015. Unveiling the species-rank abundance distribution by generalizing the Good-Turing sample coverage theory. Ecology 96:1189-1201.
Chao A. 1987. Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43:783-791.
R code for JADE: http://chao.stat.nthu.edu.tw/wordpress/paper/107_Rcode.txt
}
\author{
Youhua Chen & Tsung-Jen Shen
}
\note{
This function is a part of the original R code JADE by Chao et al. (2015) and is slightly modified for the output format.
}
\seealso{
\code{\link{DetAbu}}
}
\examples{
## As an example, Canadian-mite data are used here.	
data(CanadaMite)
## two columns represent two samples of incidence counts
X.merge = CanadaMite
## the first column is treated as the original sample
X.col1 = X.merge[,1]
Xi = X.col1
## the number of quadrats in the first sample
nT = 16
DetInc(y = Xi, nT = nT)	
}

