\name{bpec.mcmc}	
\encoding{ISO-8859-1}
\Rdversion{1.1}		
\alias{bpec.mcmc}	
\title{Markov chain Monte Carlo sampler for BPEC}	
\description{Markov chain Monte Carlo for Bayesian Phylogeographic and Ecological Clustering, implemented in C. Given a dataset of DNA sequences (non-recombinant, typically mtDNA) and their respective geographical locations (longitude and latitude), the algorithm simultaneously draws inferences about the genealogy (in the form of a haplotype tree) and the population clustering (with an unknown number of clusters). In addition, the algorithm identifies locations with high posterior probability of being ancestral. In case where additional covariates are available (e.g. climate data), these may be added to the 2-dimensional data and inserted into the analysis. }
\usage{bpec.mcmc(rawSeqs,coordsLocs,maxMig,iter,ds,postSamples,dims=-1)}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{rawSeqs}{The input DNA sequences.}	
  \item{coordsLocs}{A matrix where each row shows latitude, longitude (plus any additional covariates), plus all the haplotype numbers found at each location.  }

\item{maxMig}{The maximum number of migration events (this means that the maximum number of clusters will be \code{maxmig+1}). The number you enter here is just an upper bound, so start with \code{maxmig=6} and only increase it if you are really getting 7 clusters in return (which could mean that more clusters are appropriate). If, say, 4 clusters are needed, whether you use \code{maxMig=6} or \code{maxMig=10} (or similar), the number of clusters will collapse down to 4.  }
  \item{iter}{The number of iterations for the MCMC sampler, must be a multiple of 10. You will need quite a large number here, like 100,000. Two MCMC chains will run, after which convergence is checked. If convergence has not been reached, the output will say "NO CONVERGENCE" and you should increase the number of iterations. }
\item{ds}{This represents the parsimony relaxation parameter, with 0 being the minimum. Generally the higher \code{ds}, the more candidate trees are considered, but this comes at a computational cost. Start with \code{ds=0} and increase to \code{ds=1}, etc, observing any changes. }
\item{postSamples}{How many posterior samples per chain to save for use post-processing. A value of at least \code{PostSamples=1000} would provide a reasonable assessment of posterior uncertainty. \code{PostSamples} must not be greater than \code{iter/10}. Also, only up to 20 saved (thinned) samples are used in the \code{bpec.contourPlot} function. }
\item{dims}{The dimension, 2 for purely geographical data, +1 for each covariate (for example if environmental or phenotypic characteristics are also available).}
}
\details{This is the main command for BPEC, the details of which can be found in the provided references. In short, the model is such that, for a given haplotype tree, clusterings are the result of the migration of an individual, such that all the descendants of that individual belong to the new cluster (unless they subsequently migrate themselves). The number of migrating individuals is itself a parameter. Conditional on the clustering allocations for each sequence (note: not all samples of a haplotype need to belong to the same cluster), the geographical (and covariate) distribution for each cluster is assumed to be Gaussian. The distributions of longitude and latitude in each cluster are assumed to have an unknown covariance, whereas additional covariates are assumed to be independent. 

Posterior samples are obtained through Markov chain Monte Carlo (MCMC) sampling. Two Markov chain Monte Carlo runs are carried out. In each, the sampler sweeps through the updates of all the parameters: the haplotype tree, the root of the tree, the number of clusters, the precise clustering (as a result of migrations) the locations and variances of the clusters. Metropolis-Hastings updates are performed for most of the parameters. In the case of the number of clusters, Metropolis-Hastings Reversible Jump updates are performed. 

The space of possible clusterings is vast and highly multi-modal. The migration of an individual to a new cluster implies that all of their descendants will belong to the new cluster. This results in a combinatorial parameter space which is challenging to explore. A number of sophisticated tricks are used in order to overcome this challenge, alternating between local and global MCMC proposal moves. The biggest challenge is to converge to the region of high posterior probability in terms of the number of clusters and the cluster allocation. As such, the first 90 percent of the iterations are discarded as burn-in and only the final 10 percent are used as potential posterior samples. 

\code{bpec.mcmc} requires 2 input files in order to run:

\code{haplotypes.nex} :
   The sequence file in NEXUS format. Sequence labels should either be integers, or contain unique integers which correspond to the labels in the \code{CoordsLocsFile.txt}. For example, '1', '1_label', 'label1_label' will all be treated as haplotype 1. NOTE: BPEC will currently ignore unknown nucleotides in the inference.  
  
\code{coordsLocsFile.txt} :
   For each observation, the coordinates (latitude and longitude, please use a +/- to  indicate W or E), any other environmental or phenotypic covariates (the latitude and longitude MUST come first), plus the ID number of the haplotype (must match the number in the sequence NEXUS file). If more than one haplotype were found at a single location, these can be entered one after the other, eg:

\preformatted{
  36.88	-5.42	 24	 25
  37.00	-3.98	245	251	243	142	143	244	246	247
  43.35  1.48	153		
}
   so, in the first location (lat/long 36.88, -5.42) you have 2 sampled individuals with   haplotypes 24,25, in the second location  eight individuals etc. Sequences don't necessarily need to be collapsed onto haplotypes, the program should take  care of it.

}
\value{\code{bpec} object which can be analysed and summarised using the accessor functions \code{input()}, \code{preproc()}, \code{output.clust()}, \code{output.tree()}, \code{output.mcmc()}.}

\note{\code{bpec.mcmc} uses  \code{cexcept.h 2.0.1} (an interface for exception-handling in ANSI C) developed by Adam M. Costello and Cosmin Truta.
}

\references{
I. Manolopoulou and B.C. Emerson (2012). Phylogeographic ancestral inference using the coalescent model on haplotype trees. \emph{Journal of Computational Biology}, 19(6), 745-755.

I. Manolopoulou, L. Legarreta, B.C. Emerson, S. Brooks, and S. Tavaré (2011). A Bayesian approach to phylogeographic clustering. \emph{Interface focus}, rsfs20110054.

 Adam M. Costello and Cosmin Truta (2008) \code{cexcept.h} exception handling interface in C http://www.nicemice.net/cexcept/.

}
\author{Ioanna Manolopoulou}
\examples{
#to use the example dataset:
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
coordsLocs <- MacrocnemisCoordsLocs
rawSeqs <- MacrocnemisRawSeqs

##to use your own dataset 
#rawSeqs <- bpec.loadSeq('Haplotypes.nex')
#coordsLocs <- bpec.loadCoords("coordsLocsFile.txt")

## to set phenotypic/environmental covariate names manually, use (as appropriate)
# colnames(coordsLocs)[1:dims] <- c('lat','long','cov1','cov2','cov3')   
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 

#to run the MCMC sampler: 
bpecout <- bpec.mcmc(rawSeqs, coordsLocs, maxMig = 2, iter = 50, ds = 0, postSamples = 5, dims = 8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{phylogeography} 
\keyword{clustering}
\keyword{Bayesian}
\keyword{MCMC}
\keyword{convergence}



