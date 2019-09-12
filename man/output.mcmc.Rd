\name{output.mcmc}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{output.mcmc.bpec}
\alias{output.mcmc}
\title{Load location data}
\description{Provides various MCMC tuning parameters, as well as posterior samples for convergence assessment. }
\usage{output.mcmc(bpecout)

\method{output.mcmc}{bpec}(bpecout)
}

\arguments{
\item{bpecout}{\R object from \code{bpec.mcmc} run}		  
 }

\value{
\item{MCMCparams}{Various MCMC tuning parameters, useful for development. }
\item{codaInput}{Posterior samples from the two MCMC chains for the cluster means, cluster covariance entries, as well as the root haplotype. Note that, since the number of clusters varies from iteration to iteration, some samples are simply draws from the prior (corresponding to empty clusters). This variable can be loaded directly into the \code{coda} package for convergence analysis.  } 
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
output.mcmc(bpecout)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
