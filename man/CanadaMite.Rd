\name{CanadaMite}
\alias{CanadaMite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mite incidence in moss patches of 32 locations of western Canada (Chen et al. 2015)
}
\description{
The mite data were collected by Youhua Chen (Chen et al. 2015) in western coasts of Canada.
In Chen et al. (2015), 16 moss sampling locations were surveyed from the early days to the midst of June 2011,
while another 16 moss sampling units were surveyed from the midst of June to the early days of July 2011.

Therefore, for the dataset, it has two columns, the first column contained the incidence or occurrence information
of mites collected by first days of sampling (the early days to the midst of June 2011),
while the second column contained the incidence of mites that were collected by the last days of the sampling
}
\usage{
data(CanadaMite)
}

\references{
Chen Y, Amundrud SL, Srivastava DS (2015) Spatial variance in soil microarthropod communities: Niche, neutrality, or stochasticity? Ecoscience 21:1-14.

Chen Y, Shen TJ (2019) How many rare species in a total in additional ecological samples? In preparation.

Shen TJ, Chen YH (2018) A Bayesian weighted approach to predicting the number of newly discovered rare species. Conservation Biology, In press.

}
\author{
Youhua Chen & Tsung-Jen Shen
}

\seealso{
\code{\link{HerpetologicalData}}
}
\examples{	
data(CanadaMite)
## two columns represent two samples of species incidence data
X.merge = CanadaMite
head(X.merge)
}

