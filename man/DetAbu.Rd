\name{DetAbu}
\alias{DetAbu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Abundance-based data: the estimation of parameters for obtaining the estimation of relative abundances of observed species 
}
\description{
The estimation of parameters for obtaining the estimation of relative abundances of observed species 
}
\usage{
DetAbu(x, zero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A vector of species abundance data, i.e., the number of individuals of species 1, the number of individuals of species 2, and so forth. 
}
  \item{zero}{
A logical value for whether reserving zero frequency or not.
}
}
\value{
A vector of 3 elements: the first two values are the estimates of two parameters in Chao et al. (2015) for jointly estimating relative abundances of observed species and the third one is the estimated number of unseen species in the sample by Chao 1 estimator (Chao, 1984).
}
\references{
Chao A, Hsieh T, Chazdon R, Colwell R, Gotelli N. 2015. Unveiling the species-rank abundance distribution by generalizing the Good-Turing sample coverage theory. Ecology 96:1189-1201.

Chao A. 1984. Non-parametric estimation of the number of classes in a population. Scandinavian Journal of Statistics 11:265-270.

R code for JADE: http://chao.stat.nthu.edu.tw/wordpress/paper/107_Rcode.txt
}
\author{
Youhua Chen & Tsung-Jen Shen
}
\note{
This function is a part of the original R code JADE by Chao et al. (2015) and is slightly modified for the output format.
}
\seealso{
\code{\link{DetInc}}
}
\examples{
## As an example, Herpetological assemblage data are used here.		
data(HerpetologicalData)
## two columns represent two samples of species abundance data
X.merge = HerpetologicalData
## the first column is treated as the original sample
Xi = X.merge[,1]
DetAbu(x=Xi)	
}

