\name{bpec.treePlot}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{bpec.treePlot}
\title{Treeplot indicating cluster membership}
\description{Given a set of haplotype nodes in a tree structure, creates
  a haplotype tree plot, starting from the root at the top. The size of each haplotype node represents the
number of times it was sampled, and color represents population
cluster. Unsampled haplotypes are represented by black dots. Thickness of tree edges represents posterior certainty, with thinner edges corresponding to higher uncertainty. }
\usage{bpec.treePlot(bpecout,colorCode=c(7,5,6,3,2,8,4,9))}

\arguments{
  \item{bpecout}{\R object from \code{bpec.mcmc} run. }
   \item{colorCode}{A vector of color codes to use, ideally the same ones used in bpec.ContourPlot.}
 }

\value{

\item{graphEdgesSub}{Set of nodes and branch that can be loaded into \code{plot.igraph()}.}
\item{graphEdgesTree}{Phylogenetic tree representation that can be plotted directly using \code{plot()}.}
}

\details{The algorithm will try to plot the tree such that branches
  don't cross. The label which will be shown correspond to the original labels provided by the user. In cases where 2 sequences (with different labels) were collapsed to the same haplotype, the smallest (in number) of the two labels will be shown. }

 \references{
 I. Manolopoulou, A. Hille, B.C.Emerson B (2020). BPEC: An R Package for
Bayesian Phylogeographic and Ecological Clustering. \emph{Journal of
Statistical Software}, 92(3), 1-32. doi: 10.18637/jss.v092.i03

G. Csardi, T. Nepusz (2006). The igraph software package for complex network research. \emph{InterJournal, Complex Systems}, 1695(5), 1-9.

I. Manolopoulou and B.C. Emerson (2012). Phylogeographic ancestral inference using the coalescent model on haplotype trees. \emph{Journal of Computational Biology}, 19(6), 745-755.

G. Valiente (2009). \emph{Combinational Pattern Matching Algorithms in Computational Biology Using Perl and R}. CRC Press.
}

\author{Ioanna Manolopoulou & Axel Hille}
\examples{
#to use example dataset:
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
coordsLocs <- MacrocnemisCoordsLocs
rawSeqs <- MacrocnemisRawSeqs

##to use your own dataset 
#rawSeqs <- bpec.loadSeq('haplotypes.nex')
#coordsLocs <- bpec.loadCoords("coordsLocsFile.txt")

## to set phenotypic/environmental covariate names manually, use (as appropriate)
# colnames(coordsLocs)[1:dims] <- c('lat','long','cov1','cov2','cov3')   
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 

#to run the MCMC sampler: 

bpecout <- bpec.mcmc(rawSeqs, coordsLocs, maxMig = 2, iter = 50, ds = 0, postSamples = 1, dims = 8)
bpec.Tree <- bpec.treePlot(bpecout)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree} 
\keyword{phylogenetics}
\keyword{haplotype}
