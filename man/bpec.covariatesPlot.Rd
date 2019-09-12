\name{bpec.covariatesPlot}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{bpec.covariatesPlot}
\title{Plot posterior marginal distribution of each covariate for the inferred clusters}
\description{For each covariate, the posterior marginal distribution for each cluster is shown, with colours representing cluster index. The posterior median as well as .10 and .90 quantiles are shown on the same plot in order to provide a measure of uncertainty.  }
\usage{bpec.covariatesPlot(bpecout, colorCode=c(7,5,6,3,2,8,4,9))}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{bpecout}{\R object from \code{bpec.mcmc} run }
  \item{colorCode}{A vector of color codes}
}
\details{Unlike the contour plot, where all posterior samples are shown on the map, in this case quantiles are used pointwise along the axis of each covariate. The plots are shown in the original scale of the covariates and axis labels can be fed into the function to aid interpretation. }
\value{
Produces plots for all additional (environmental or phenotypic) covariate, where posterior medians as well as credible intervals are shown for each cluster. 
}
\author{Ioanna Manolopoulou}
\references{
I. Manolopoulou, L. Legarreta, B.C. Emerson, S. Brooks, and S. Tavaré (2011). A Bayesian approach to phylogeographic clustering. \emph{Interface focus}, rsfs20110054.
}

\examples{
#to use example dataset:
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
rawSeqs <- MacrocnemisRawSeqs
coordsLocs <- MacrocnemisCoordsLocs

##to use your own dataset 
#rawSeqs <- bpec.loadSeq('Haplotypes.nex')
#coordsLocs <- bpec.loadCoords("coordsLocsFile.txt")

## to set phenotypic/environmental covariate names manually, use (as appropriate)
# colnames(CoordsLocs)[1:dims] <- c('lat','long','cov1','cov2','cov3')   
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 

#for the analysis: 
#check the helpfile of bpec.mcmc using ?bpec.mcmc

bpecout <- bpec.mcmc(rawSeqs, coordsLocs, maxMig = 2, iter = 100, ds = 0, postSamples = 5, dims = 8)

#if there are also environmental covariates available: 
par(mfrow=c(2,3)) #this splits the plot window into 2x3 to fit all the covariates
bpec.covariatesPlot(bpecout)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{clustering}
\keyword{Bayesian}
\keyword{covariates}
\keyword{environmental}
\keyword{phenotypic}

