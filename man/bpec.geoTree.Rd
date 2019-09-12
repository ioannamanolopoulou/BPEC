\name{bpec.geoTree}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{bpec.geoTree}
\title{Phylogeographic GoogleEarth plot}
\description{Given the output of an MCMC run, outputs phylogeographic code to show migration on a map through GoogleEarth.  }
\usage{bpec.geoTree(bpecout,file="GoogleEarthTree.kml")}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bpecout}{\R object from \code{bpec.mcmc} run }
 
    \item{file}{Filename for the .kml file to be created }
}
\details{ Additionally to the output variables, a .kml file is created which can be loaded into GoogleEarth.  }

\value{
\item{phy}{Set of nodes and branches that can be loaded into \code{bpec.phylo2GE()}}
\item{geo}{Set of haplotypes and locations that can be loaded into \code{bpec.phylo2GE()}}
}

\references{
I. Manolopoulou and B.C. Emerson (2012). Phylogeographic ancestral inference using the coalescent model on haplotype trees. \emph{Journal of Computational Biology}, 19(6), 745-755.

G. Valiente (2009). \emph{Combinational Pattern Matching Algorithms in Computational Biology Using Perl and R}. CRC Press.

N. Arrigo, L.P. Albert, P.G. Mickelson and M.S. Barker (2012). Quantitative visualization of biological data in Google Earth using R2G2, an R CRAN package. \emph{Molecular Ecology Resources}, 12(6), 1177-1179. 

}

\note{\code{bpec.phylo2GE} uses the function \code{curvy.R} from the archived package \code{R2G2} of Arrigo et al (2012). 
}


\author{Axel Hille}
\examples{
\dontrun{
#to use example dataset:
data(MacrocnemisRawSeqs)
data(MacrocnemisCoordsLocs)
coordsLocs <- MacrocnemisCoordsLocs
rawSeqs <- MacrocnemisRawSeqs

## to use your own dataset
# rawSeqs <- bpec.loadSeq('Haplotypes.nex')
# coordsLocs <- bpec.loadCoords("coordsLocsFile.txt")

## to set phenotypic/environmental covariate names manually, use (as appropriate)
# colnames(coordsLocs)[1:dims] <- c('lat','long','cov1','cov2','cov3')   
## where dims is the corresponding number of measurements available 
## (2 for latitude and longitude only, add one for each additional available measurement) 

#to run the MCMC sampler: 
bpecout <- bpec.mcmc(rawSeqs, coordsLocs, maxMig = 2, iter = 50, ds = 0, postSamples = 5, dims = 8)
bpec.Geo <- bpec.geoTree(bpecout, file = "GoogleEarthTree.kml")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{tree} 
\keyword{phylogenetics}
\keyword{haplotype}

