\name{bpec.contourPlot}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{bpec.contourPlot}
\title{Plot contours of inferred clusters on geographical map}
\description{Plots the contours for a set of population clusters, adding
points for sampling locations, with larger points representing inferred
ancestral locations. Assuming that the data input is in latitude and
longitude form, the contour is superimposed onto the map of the
corresponding region.    }
\usage{bpec.contourPlot(bpecout, GoogleEarth = 0, 
colorCode = c(7,5,6,3,2,8,4,9), mapType = 'plain', mapCentre = NULL, zoom = 6)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bpecout}{\R object from \code{bpec.mcmc} run }
\item{GoogleEarth}{If 1, .kml files are produced which can be opened with GoogleEarth.}
 \item{colorCode}{(optional) A vector of color codes.}
 \item{mapType}{Options are 'none' (no map shown in the background), 'plain' (outline shown only), 'google' (google maps), 'osm' (openstreetview maps).}
\item{mapCentre}{(optional) The longitude and latitude to use as the centre of the map. If not provided, the midpoint of sampling range will be used.  }
\item{zoom}{The level of zoom into the map, default 6.}
 }
\details{Contours with level 0.5 of the posterior mean centres amd covariances of the clusters are shown on the map, with colour representing cluster index. 
Posterior samples are shown with semi-transparent filled contours, so that uncertainty can be represented directly on the map together with the posterior means. The sampling locations are shown as black dots, with the top three (in terms of posterior probability) ancestral locations represented by larger dots. 
If \code{GoogleEarth=1}, this same plot is exported in files that can be directly loaded into Google Earth. 
 }
\value{
Produces a contour plot overlayed on a map. If \code{GoogleEarth=1}, files which may be loaded directly into Google Earth are produced. 
}
\author{Ioanna Manolopoulou, Axel Hille and Steve Brooks}
\references{
I. Manolopoulou, L. Legarreta, B.C. Emerson, S. Brooks, and S. Tavaré (2011). A Bayesian approach to phylogeographic clustering. \emph{Interface focus}, rsfs20110054.

S.P. Brooks, I. Manolopoulou, and B.C. Emerson (2007). \emph{Assessing the Effect of Genetic Mutation - A Bayesian Framework for Determining Population History from DNA Sequence Data}. Bayesian Statistics 8. Oxford University Press. 
}

\examples{
#to use example dataset:
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
rawSeqs <- MacrocnemisRawSeqs
coordsLocs <- MacrocnemisCoordsLocs

##to use your own dataset: 
#rawSeqs <- bpec.loadSeq('Haplotypes.nex')
#coordsLocs <- bpec.loadCoords("coordsLocsFile.txt")

## to set phenotypic/environmental covariate names manually, use (as appropriate)
# colnames(coordsLocs)[1:dims] <- c('lat','long','cov1','cov2','cov3')   
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 

#for the analysis: 
#check the helpfile of bpec.mcmc using ?bpec.mcmc
colorCode <- c(7,5,6,3,2,8) #default colour scheme
bpecout <- bpec.mcmc(rawSeqs, coordsLocs, maxMig = 2, iter = 20, ds = 0, postSamples = 2, dims = 8)
bpec.contourPlot(bpecout, GoogleEarth = 0, colorCode, mapType = 'plain') 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{map} 
\keyword{Bayesian}
\keyword{contours}

