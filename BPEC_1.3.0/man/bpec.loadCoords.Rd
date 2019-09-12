\name{bpec.loadCoords}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{bpec.loadCoords}
\title{Load location data}
\description{Create \R objects of  locations from  corresponding file. }
\usage{bpec.loadCoords(coordsFile, header = FALSE)}

\arguments{
\item{coordsFile}{The name of the coordinate and sequence file in full, eg "coordsLocs.txt".}
\item{header}{If TRUE, then the first row of the file is assumed to be variable names. }
 }

\value{
\item{coordsLocs}{\R object of locations and corresponding haplotypes.}
}

\details{Given the names of the file, this command will create an \R objects of the coordinates and locations that can be fed into \code{bpec.mcmc}. }


\author{Ioanna Manolopoulou & Axel Hille}
\examples{

#to load existing example file from BPEC
coordsFile <- system.file("coordsLocsFile.txt",package = "BPEC")
coordsLocs <- bpec.loadCoords(coordsFile, header = TRUE)

\dontrun{
#to load a different input file "coordsLocs.txt"
coordsLocs <- bpec.loadCoords("coordsLocs.txt", header = TRUE)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
