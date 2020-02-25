\name{MacrocnemisCoordsLocsMini}
\encoding{ISO-8859-1}
     \docType{data}
     \alias{MacrocnemisCoordsLocsMini}
     \title{Coordinates and haplotypes for the R.macrocnemis Brown Frog dataset}
     \description{A reduced version of the 
\emph{Rana macrocnemis} (Amphibia, Ranidae) dataset, the Near Eastern brown frog, from a transition zone of landscapes in Georgia, Caucasus.
       Contains the coordinates and ecological characteristics of each sampling location along with the haplotypes found in each location. It is called using \code{data(MacrocnemisCoordsLocsMini)}.
     }
     \usage{MacrocnemisCoordsLocsMini}
     \format{In total 21 sampling locations with one haplotype in each. For each observation, longitude, latitude, plus Annual Mean 
Temperature (degrees Celsius x 10) from the \url{http://www.worldclim.org} database.}
    
     \references{
D. Tarkhnishvili, A. Hille  and W. B?hme (2001) 
Humid forest refugia, speciation and secondary introgression between evolutionary lineages: differentiation in a Near Eastern brown frog, Rana macrocnemis. Biological Journal of the Linnean Society.

      A. Hille, I. Manolopoulou and D. Tarkhnisvili (2015)
Detecting landscape-dependent selection with Bayesian phylogeographic and ecological clustering: 
comparison of recombinant microsatellite and mtDNA haplotype data (unpublished).
     }
     \examples{
#to use the example dataset:
data(MacrocnemisCoordsLocsMini)
coordsLocs <- MacrocnemisCoordsLocsMini
}
     \keyword{Brown Frog}
