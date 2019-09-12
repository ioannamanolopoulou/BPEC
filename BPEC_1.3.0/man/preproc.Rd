\name{preproc}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{preproc.bpec}
\alias{preproc}
\title{Load location data}
\description{Provides various MCMC tuning parameters, as well as posterior samples for convergence assessment. }
\usage{preproc(bpecout)

\method{preproc}{bpec}(bpecout)
}

\arguments{
\item{bpecout}{\R object from \code{bpec.mcmc} run}
		  
 }

\value{
\item{seq}{The output DNA sequences of distinct haplotypes, collapsed to effective nucleotide sites (both sampled and missing sequences which were inferred).}
\item{seqsFile}{A vector of the numerical labels of each haplotype.}
\item{seqLabels}{Correspondence vector for each of the processed observations to the original haplotype labels. }
\item{seqIndices}{Correspondence vector for each of the original observations to the resulting haplotype labels.}
\item{seqLength}{The effective length of the input sequences, given by the number of variable nucleotide sites which are informative. In other words, if two different nucleotide sites are variable in exactly the same haplotypes, then they effectively provide information of a single site.}
\item{noSamples}{The number of times each haplotype was observed in the sample. }
\item{count}{The number of output sequences.}
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
preproc(bpecout)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
