\name{MacrocnemisCoordsLocs}
\encoding{ISO-8859-1}
     \docType{data}
     \alias{MacrocnemisCoordsLocs}
     \title{Coordinates and haplotypes for the R.macrocnemis Brown Frog dataset}
     \description{
\emph{Rana macrocnemis} (Amphibia, Ranidae) dataset, the Near Eastern brown frog, from a transition zone of landscapes in Georgia, Caucasus.
       Contains the coordinates and ecological characteristics of each sampling location along with the haplotypes found in each location. It is called using \code{data(MacrocnemisCoordsLocs)}.
     }
     \usage{MacrocnemisCoordsLocs}
     \format{In total 40 sampling locations with one haplotype in each. For each observation, longitude, latitude, plus 6 additional environmental coordinates are available: Four bioclimatic 
variables from the \url{http://www.worldclim.org} database (Annual Mean 
Temperature (degrees Celsius x 10), Temperature Annual 
Range (100 x standard deviation of monthly mean temperature), 
Annual Precipitation (in mm),  Precipitation 
Seasonality), Coefficient of Variation (CV), altitude in meters 
as a proxy for a digital elevation model (DEM) and the land cover 
map (GLC2000) from \url{http://worldgrids.org}. In fact some sampling locations may be identical in terms of longitude and latitude, but are differentiated because of their additional covariates. }
    
     \references{
D. Tarkhnishvili, A. Hille  and W. B?hme (2001) 
Humid forest refugia, speciation and secondary introgression between evolutionary lineages: differentiation in a Near Eastern brown frog, Rana macrocnemis. Biological Journal of the Linnean Society.

      A. Hille, I. Manolopoulou and D. Tarkhnisvili (2015)
Detecting landscape-dependent selection with Bayesian phylogeographic and ecological clustering: 
comparison of recombinant microsatellite and mtDNA haplotype data (unpublished).
     }
     
       \examples{
#to use the example dataset:
data(MacrocnemisCoordsLocs)
coordsLocs <- MacrocnemisCoordsLocs
}
     \keyword{Brown Frog}
