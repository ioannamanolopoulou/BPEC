\name{output.clust}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{output.clust.bpec}
\alias{output.clust}
\title{Load location data}
\description{Posterior output for the clustering parameters. }
\usage{output.clust(bpecout)

\method{output.clust}{bpec}(bpecout)
}

\arguments{
\item{bpecout}{\R object from \code{bpec.mcmc} run}
		  
 }

\value{
\item{sampleMeans}{A set of posterior samples of the cluster means (i.e. centres).}
\item{sampleCovs}{A set of posterior samples of the cluster covariances (i.e. shapes).}
\item{sampleIndices}{A set of posterior samples of the cluster allocations of each observation.} 
\item{clusterProbs}{For each haplotype, posterior probabilities that it belongs to each cluster.}
}

\author{Ioanna Manolopoulou & Axel Hille}
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
output.clust(bpecout)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
