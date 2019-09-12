\name{input}
\encoding{ISO-8859-1}
\Rdversion{1.2}
\alias{input.bpec}
\alias{input}
\title{Load location data}
\description{Shows all input files and settings used in a BPEC run. }
\usage{input(bpecout)

\method{input}{bpec}(bpecout)
}

\arguments{
\item{bpecout}{\R object from \code{bpec.mcmc} run}		  
 }

\value{
\item{seqCountOrig}{The number of input sequences.}
\item{seqLengthOrig}{The length of the input sequences.} 
\item{iter}{The number of MCMC iterations. }
\item{ds}{The parsimony relaxation parameter.} 
\item{coordsLocs}{The input coordinate and observation file. }
\item{coordsDims}{The input dimension (2 for purely geographical data).}
\item{locNo}{The number of distinct sampling locations. }
\item{locData}{The list of coordinates of each observation.}  
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
input(bpecout)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
