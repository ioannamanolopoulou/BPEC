\name{bpec}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{bpec}
\alias{plot.bpec}
\alias{summary.bpec}
\alias{mean.bpec}
\alias{print.bpec}

\title{Bayesian Phylogeographic and Ecological Clustering (BPEC)}
		
\description{	
Bayesian Phylogeographic and Ecological Clustering (BPEC) is aimed at identifying genetically, geographically and ecologically distinct population clusters and drawing inferences about ancestral locations. Given a dataset of DNA sequences (non-recombinant, typically mtDNA) and their respective geographical locations (longitude and latitude) along with additional environmental or phenotypic characteristics, the algorithm simultaneously draws inferences about the genealogy (in the form of a haplotype tree) and the population clustering. In addition, the algorithm identifies locations with high posterior probability of being ancestral. 
}

\usage{
bpec(seqsFile, coordsFile, dims = 2, iter = 100000, postSamples = 100, 
     maxMig = 5, ds = 0, colorCode = c(7,5,6,3,2,8))

\method{plot}{bpec}(x, GoogleEarth = 0, colorCode = c(7,5,6,3,2,8),...)
\method{summary}{bpec}(object,...)
\method{print}{bpec}(x,...)
\method{mean}{bpec}(x,...)
}


\arguments{
   \item{x}{\R object from \code{bpec.mcmc} run. }
   \item{object}{\R object from \code{bpec.mcmc} run. }
  
  \item{colorCode}{A vector of color codes to use, ideally the same ones used in bpec.contourPlot.}
\item{seqsFile}{The name of the NEXUS file in full, eg "SeqsFile.nex".}
\item{coordsFile}{the name of the coordinate and sequence file in full, eg "CoordsLocs.txt".}
\item{maxMig}{The maximum number of migration events (this means that the maximum number of clusters will be \code{maxmig+1}). The number you enter here is just an upper bound, so start with \code{maxmig=6} and only increase it if you are really getting 7 clusters in return (which could mean that more clusters are appropriate). If, say, 4 clusters are needed, whether you use \code{maxmig=6} or \code{maxmig=10} (or similar), the number of clusters will collapse down to 4.  }
  \item{iter}{The number of iterations for the MCMC sampler, must be a multiple of 10. You will need quite a large number here, like 100,000. Two MCMC chains will run, after which convergence is checked. If convergence has not been reached, the output will say "NO CONVERGENCE" and you should increase the number of iterations. }
\item{ds}{This represents the parsimony relaxation parameter, with 0 being the minimum. Generally the higher \code{ds}, the more candidate trees are considered, but this comes at a computational cost. Start with \code{ds=0} and increase to \code{ds=1}, etc, observing any changes. }
\item{postSamples}{How many posterior samples per chain to save for use post-processing. A value of at least \code{PostSamples=1000} would provide a reasonable assessment of posterior uncertainty. \code{PostSamples} must not be greater than \code{iter/10}. Also, only up to 20 saved (thinned) samples are used in the \code{bpec.contourPlot} function. }
\item{dims}{The dimension, 2 for purely geographical data, +1 for each covariate (for example if environmental or phenotypic characteristics are also available).}
\item{GoogleEarth}{If 1, .kml files are produced which can be opened with GoogleEarth.}
\item{...}{Other default options.}
}

\value{ \code{bpec} object which can be analysed using \code{input()}, \code{data()}, \code{output.clust()}, \code{output.tree()} and \code{output.mcmc()}.}


\details{
bpec requires 2 input files in order to run:

\code{haplotypes.nex} :
  The sequence file in NEXUS format. Sequence labels should either be integers, or contain unique integers which correspond to the labels in the \code{CoordsLocsFile.txt}. For example, '1', '1_label', 'label1_label' will all be treated as haplotype 1. NOTE: bpec will currently ignore unknown nucleotides in the inference.  

\code{coordsLocsFile.txt} :
   For each observation, the coordinates (latitude and longitude, please use a +/- to  indicate W or E), any other environmental or phenotypic covariates (the latitude and longitude MUST come first), plus the ID number of the haplotype (must match the number in the sequence NEXUS file). If more than one haplotype were found at a single location, these can be entered one after the other, eg:

\preformatted{
  36.88	-5.42	 24	 25
  37.00	-3.98	245	251	243	142	143	244	246	247
  43.35  1.48	153		
}
   so, in the first location (lat/long 36.88, -5.42) you have 2 sampled individuals with   haplotypes 24,25, in the second location  eight individuals etc. Sequences don't necessarily need to be collapsed onto haplotypes, the program should take  care of it.
 
 
The main function is \code{bpec.mcmc} and runs a Markov chain Monte Carlo sampler in order to obtain posterior samples for all the parameters simultaneously: the haplotype tree, the ancestral nodes, the number of population clusters as well as their means and covariances. 
 
 
Three plotting functions are available. \code{bpec.contourPlot} shows the inferred population clusters superimposed on a world map. \code{bpec.treePlot} shows the Maximum A Posteriori rooted haplotype tree, indicating posterior cluster membership and number of times each haplotype was sampled.  It can also output .kml files which can be loaded into Google Earth.  \code{bpec.CovariatesPlot} shows the posterior density estimates for each additional covariate of each cluster.  
}

\note{
\code{bpec.mcmc} uses  \code{cexcept.h 2.0.1} (an interface for exception-handling in ANSI C) developed by Adam M. Costello and Cosmin Truta.
}

\author{
Ioanna Manolopoulou <ioanna.manolopoulou@gmail.com>, Axel Hille <axel.hille@gmx.net>
}
\references{
I. Manolopoulou, A. Hille, B.C.Emerson B (2020). BPEC: An R Package for
Bayesian Phylogeographic and Ecological Clustering. \emph{Journal of
Statistical Software}, 92(3), 1-32. doi: 10.18637/jss.v092.i03

I. Manolopoulou and B.C. Emerson (2012). Phylogeographic ancestral inference using the coalescent model on haplotype trees. \emph{Journal of Computational Biology}, 19(6), 745-755.

I. Manolopoulou, L. Legarreta, B.C. Emerson, S. Brooks, and S. Tavar? (2011). A Bayesian approach to phylogeographic clustering. \emph{Interface focus}, rsfs20110054.

S.P. Brooks, I. Manolopoulou, and B.C. Emerson (2007). \emph{Assessing the Effect of Genetic Mutation - A Bayesian Framework for Determining Population History from DNA Sequence Data}. Bayesian Statistics 8. Oxford University Press. 

 Adam M. Costello and Cosmin Truta (2008) \code{cexcept.h} exception handling interface in C, available at website http://www.nicemice.net/cexcept/.

}

\keyword{phylogeography} 
\keyword{clustering}
\keyword{Bayesian}
\keyword{phylogenetics}

\examples{
  ## if you want to load the `mini' example Brown Frog dataset
  data(MacrocnemisRawSeqs)
  data(MacrocnemisCoordsLocsMini)
  rawSeqs <- MacrocnemisRawSeqs
  coordsLocs <- MacrocnemisCoordsLocsMini

  dims <- 3 #this is 2 if you only have geographical longitude/latitude. 
#(add 1 for each environmental or phenotypic covariate)
maxMig <- 2 #you will need a higher maximum number of migrations, suggest 7
ds <- 0 #start with ds=0 and increase to 1 and then to 2
iter <- 1000 #you will need far more iterations for convergence, start with 100,000
postSamples <- 100 #you will need at least 100 saved posterior samples

#run the Markov chain Monte Carlo sampler
bpecout <- bpec.mcmc(rawSeqs,coordsLocs,maxMig,iter,ds,postSamples,dims)

par(mar=c(0,0,0,0),pty="m",mfrow=c(1,2))  #no plot margins, plot contours and tree side-by-side
# plot geographical cluster contour map
bpec.contourPlot(bpecout,GoogleEarth=0) 

# plot tree network with cluster indicators
bpec.Tree <- bpec.treePlot(bpecout)   

# now also plot the environmental covariates
 bpec.covariatesPlot(bpecout) 
 
 bpec.Geo <- bpec.geoTree(bpecout,file="GoogleEarthTree.kml")


\dontrun{
# if you want to load the example burnet moth dataset
data(TransalpinaRawSeqs)
data(TransalpinaCoordsLocs)
rawSeqs <- TransalpinaRawSeqs
coordsLocs <- TransalpinaCoordsLocs

##if you want to use your own dataset, use setwd() to enter the correct folder, 
##then run the command below, changing the input parameters if necessary
 #rawSeqs <- bpec.loadSeq('haplotypes.nex')
 #coordsLocs <- bpec.loadCoords("coordsLocsFile.txt")

 ## to set phenotypic/environmental covariate names manually, use (as appropriate)
# colnames(CoordsLocs)[1:dims] <- c('lat','long','cov1','cov2','cov3')  
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 
 
dims <- 2 #this is 2 if you only have geographical longitude/latitude. 
#(add 1 for each environmental or phenotypic covariate)
maxMig <- 5 #you will need a higher maximum number of migrations, suggest 7
ds <- 0 #start with ds=0 and increase to 1 and then to 2
iter <- 10000 #you will need far more iterations for convergence, start with 100,000
postSamples <- 2 #you will need at least 100 saved posterior samples

#run the Markov chain Monte Carlo sampler
bpecout <- bpec.mcmc(rawSeqs,coordsLocs,maxMig,iter,ds,postSamples,dims)

par(mar=c(0,0,0,0),pty="m",mfrow=c(1,2)) #No plot margins. Contours and tree side-by-side
# plot geographical cluster contour map
bpec.contourPlot(bpecout, GoogleEarth=0, mapType = 'plain') 

# plot tree network with cluster indicators
bpec.Tree <- bpec.treePlot(bpecout)  

## if you want to load the example Brown Frog dataset
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
rawSeqs <- MacrocnemisRawSeqs
coordsLocs <- MacrocnemisCoordsLocs

dims <- 8 #this is 2 if you only have geographical longitude/latitude. 
#(add 1 for each environmental or phenotypic covariate)
maxMig <- 4 #you will need a higher maximum number of migrations, suggest 7
ds <- 2 #start with ds=0 and increase to 1 and then to 2
iter <- 10000 #you will need far more iterations for convergence, start with 100,000
postSamples <- 2 #you will need at least 100 saved posterior samples

#run the Markov chain Monte Carlo sampler
bpecout <- bpec.mcmc(rawSeqs,coordsLocs,maxMig,iter,ds,postSamples,dims)

par(mar=c(0,0,0,0),pty="m",mfrow=c(1,2))  #no plot margins, plot contours and tree side-by-side
# plot geographical cluster contour map
bpec.contourPlot(bpecout,GoogleEarth=0) 

# plot tree network with cluster indicators
bpec.Tree <- bpec.treePlot(bpecout)   

# now also plot the environmental covariates
 par(mfrow=c(2,3)) #split the plot window into 2x3 to fit all the covariates
 bpec.covariatesPlot(bpecout) 
 
 bpec.Geo <- bpec.geoTree(bpecout,file="GoogleEarthTree.kml")
 }



}
