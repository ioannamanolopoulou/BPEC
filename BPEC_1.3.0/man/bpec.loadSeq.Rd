\name{bpec.loadSeq}
\encoding{ISO-8859-1}
\Rdversion{1.1}
\alias{bpec.loadSeq}
\title{Load sequence data}
\description{Create \R objects of sequences from NEXUS file. }
\usage{bpec.loadSeq(seqsFile)}

\arguments{
 \item{seqsFile}{The name of the NEXUS file in full, eg "SeqsFile.nex".}
 }

\value{

\item{rawSeqs}{\R object of sequences.}
}

\details{Given the names of the two files, this command will create an \R objects of the sequences that can be fed into \code{bpec.MCMCout}. }


\author{Ioanna Manolopoulou & Axel Hille}
\examples{
#to load existing example file from BPEC
seqsFile <- system.file("haplotypes.nex",package = "BPEC")
rawSeqs <- bpec.loadSeq(seqsFile)

\dontrun{
#to load a different input file "haplotypes.nex"
rawSeqs <- bpec.loadSeq('haplotypes.nex')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
