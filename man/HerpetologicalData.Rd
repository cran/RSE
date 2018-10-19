\name{HerpetologicalData}
\alias{HerpetologicalData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Abundance of herpetofauna in the conserved and human disturbed areas of Mexico (Suazo-Ortuno et al. 2008)
}
\description{
Suazo-Ortuno et al. (2008) studied how the conversion of tropical forest to agricultural mosaic influenced herpetofaunal distribution and community structure
in conserved and human disturbed forest areas of neotropical Mexico.

Therefore, the dataset used here, the first and second columns represented species abundance 
in the conserved and disturbed areas, respectively. 
}
\usage{
data(HerpetologicalData)
}

\references{
Suazo-Ortuno I, Alvarado-Diaz J, Martines-Ramos M (2008) Effects of conversion of dry tropical forest to agricultural mosaic on herpetofaunal assemblages. Conservation Biology 22: 362-374.

Chen Y, Shen TJ (2019) How many rare species in a total in additional ecological samples? In preparation.

Shen TJ, Chen YH (2018) A Bayesian weighted approach to predicting the number of newly discovered rare species. Conservation Biology, In press.

}
\author{
Youhua Chen & Tsung-Jen Shen
}

\seealso{
\code{\link{CanadaMite}}
}
\examples{	
data(HerpetologicalData)
## two columns represent two samples of species abundance data
X.merge = HerpetologicalData
head(X.merge)
}

