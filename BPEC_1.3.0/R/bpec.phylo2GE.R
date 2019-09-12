# Citation:
# 1
# Combinational Pattern Matching Algorithms in Computational Biology
# Using Perl and R
# by Gabriel Valiente
# CRC Press, Chapman & Hall, Boca Raton 2009
# R source code: http://www.lsi.upc.edu/~valiente/comput-biol/

# 2
#Quantitative visualization of biological data in Google Earth
#using R2G2, an R CRAN package
#Nils Arrigo, Loren P. Albert, Pascal G. Mickelson and Michael S. Barker
#Molecular Ecology Resources (2012) doi: 10.1111/1755-0998.12012


# Usage: Using the treeEdges, creating a haplotype network in Newick format,
# converting the haplotype network to input file "phy" and re-arranging the
# coordsLocsFile to input "geo" for final use in the package "R2G2";
# KML produced using the Phylo2GE function.

curvy <-
function(f, startDD, stopDD){  
  radian = 180 / pi

  lat1 = startDD[2] / radian
  lon1 = startDD[1] / radian 
  lat2 = stopDD[2] / radian
  lon2 = stopDD[1] / radian

  d = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1-lon2))
  A = sin((1 - f) * d) / sin(d)
  B = sin(f * d) / sin(d)
  x = as.numeric(A * cos(lat1) * cos(lon1) +  B * cos(lat2) * cos(lon2))
  y = as.numeric(A * cos(lat1) * sin(lon1) +  B * cos(lat2) * sin(lon2))
  z = as.numeric(A * sin(lat1) +  B * sin(lat2))
  lat=atan2(z, sqrt(x^2 + y^2))
  lon=atan2(y, x)
  rev(c(lat * 180/pi, lon * 180/pi))
  }


mynei <- function(net,v,mode){
  if(mode=="out"){
    output = setdiff(unlist(neighborhood(net, 1, nodes=v, mode="out")), v)
  }
  if(mode=="in"){
    output = setdiff(unlist(neighborhood(net, 1, nodes=v, mode="in")), v)
  }
  return(output)
}

network.to.newick <- function (net) {
  #  r <- V(net)[which(degree(net,mode="in")==0)-1]
  r = V(net)[which(degree(net,mode="in")==0)]
  newick = obtain.newick(net,r,c())
  return(paste(c(newick,";"),collapse=""))
}

obtain.newick <- function (net,v,newick) {
  if (length(V(net)[mynei(net,v,"out")]) != 0) {
    if ( length(V(net)[mynei(net,v,"in")]) != 0
         && match(v,V(net)[mynei(net,V(net)[mynei(net,v,"in")],"out")]) != 1 )
      newick = c(newick,",")
    newick = c(newick,"(")
  }
  for (w in V(net)[mynei(net,v,"out")]) {
    newick = obtain.newick(net,w,newick)
  }
  if (length(V(net)[mynei(net,v,"out")]) != 0)
    newick = c(newick,")")
  if ( length(V(net)[mynei(net,v,"in")]) != 0
       && length(V(net)[mynei(net,v,"out")]) == 0
       && match(v,V(net)[mynei(net,V(net)[mynei(net,v,"in")],"out")]) != 1 )
    newick = c(newick,",")
  c(newick,V(net)[v]$name)
}

network.to.tree <- function (net) {
  V(net)$del = FALSE
  E(net)$del = FALSE
  #  r <- V(net)[which(degree(net,mode="in")==0)-1]
  r = V(net)[which(degree(net,mode="in")==0)]
  net = postorder.traversal(net,r)
  net = delete.vertices(net,V(net)[V(net)$del])
  net = delete.edges(net,E(net)[E(net)$del])
  network.to.newick(net)
}

postorder.traversal <- function (net,v) {
  for (w in V(net)[mynei(net,v,"out")]) {
    net <- postorder.traversal(net,w)
  }
  VW = E(net)["from" == v & E(net)$del==FALSE]
  if (length(c(VW)) == 1) {
    u = V(net)[mynei(net,v,"in")]
    w = get.edge(net,VW)[2]
    E(net,path=c(V(net)[u],V(net)[v]))$del = TRUE
    E(net,path=c(V(net)[v],V(net)[w]))$del = TRUE
    V(net)[v]$del = TRUE
    net = add.edges(net,c(V(net)[u],V(net)[w]))
    E(net,path=c(V(net)[u],V(net)[w]))$del = FALSE
  }
  net
}

explode <- function (net) {
  trees = c()
  #  H <- V(net)[which(degree(net,mode="in")>1)-1]
  H = V(net)[which(degree(net,mode="in")>1)]
  if (length(H) == 0) {
    trees = c(trees,network.to.tree(net))
  } else {
    v = c(H)[1]
    for (u in V(net)[mynei(net,v,"in")]) {
      new = delete.edges(net,E(net,P=c(u,v)))
      trees = c(trees,explode(new))
    }
  }
 # print(trees)
}


bpec.phylo2GE <-
function (geo, phy, resol = 0.05, minAlt = 1000, maxAlt = 2.5e+05, 
    goo = "GoogleEarthTree.kml") 
{    
    geo = geo[match(phy$tip.label, geo[, 1]), ]
 
    Ntaxa = length(phy$tip.label)
    Ntot = Ntaxa + phy$Nnode
    geo.childs = data.frame(node = 1:nrow(geo), geo[, 2:3])
    for (i in Ntot:(Ntaxa + 1)) {
        childs = phy$edge[phy$edge[, 1] == i, 2]
        desc = geo.childs[match(childs, geo.childs$node), 2:3]
        if (length(childs) > 2) {
            anc.geo = c(i, colMeans(desc))
        }
        else {
            if (sum(desc[1, ] == desc[2, ]) == 0) {
                anc.geo = c(i, curvy(0.5, desc[1, ], desc[2, 
                  ]))
            }
            else {
                anc.geo = c(i, as.numeric(desc[1, ]))
            }
        }
        geo.childs = rbind(geo.childs, anc.geo)
    }
    gtmp = geo.childs[order(geo.childs[, 1]), 2:3]
    root = phy$edge[1, 1]
    edgenum = nrow(phy$edge)
    ages = NULL
    for (i in 1:edgenum) {
        if (phy$edge[i, 1] == root) 
            tmp = 0
        else {
            tmp = ages[which(phy$edge[, 2] == phy$edge[i, 1])]
        }
        ages = c(ages, tmp + phy$edge.length[i])
    }
    tmp = rbind(c(Ntaxa + 1, 0), cbind(phy$edge[, 2], ages))
    tmp = tmp[order(tmp[, 1]), ]
    meta = data.frame(node = tmp[, 1], Lon = gtmp[, 1], Lat = gtmp[, 
        2], age = tmp[, 2])
    meta$age = 1 - (meta$age/(max(meta$age)))
    meta$age = minAlt + meta$age * maxAlt
    cat("<?xml version=\"1.0\"?>\n<kml xmlns=\"http://earth.google.com/kml/2.0\">\n<Document>", 
        file = goo, append = FALSE)
    cat("<description>Produced using bpec and Phylo2GE</description>\n<name>R2G2 - Phylo2GE</name>\n<open>0</open>", 
        file = goo, append = TRUE)
    for (i in 1:Ntaxa) {
        xyz = as.numeric(meta[i, 2:4])
        cat("<Placemark>\n<name>", phy$tip.label[i], "</name>\n<LookAt>\n<longitude>", 
            xyz[1], "</longitude>\n<latitude>", xyz[2], "</latitude>\n<range>", 
            xyz[3], "</range>\n</LookAt>", file = goo, append = TRUE)
    cat("<Style>\n<IconStyle>\n<scale>1.2</scale>\n<Icon>\n<href>http://maps.google.com/mapfiles/kml/paddle/red-circle.png</href>\n</Icon>\n</IconStyle>\n</Style>\n
        <Point>\n<altitudeMode>relativeToGround</altitudeMode>\n<extrude>1</extrude>\n<coordinates>", 
            paste(xyz, collapse = ",", sep = " "), "</coordinates>\n</Point>\n</Placemark>\n", 
            file = goo, append = TRUE)
    }
    cat("<Style id=\"unselectedLine\">\n<LineStyle>\n<color>ff2fffad</color>\n<width>3</width>\n</LineStyle>\n</Style>\n<Style id=\"selectedLine\">\n<LineStyle>\n<color>ff2fffad</color>\n<width>4</width>\n</LineStyle>\n</Style><Folder><name>Edges</name>", 
        file = goo, append = TRUE)
    for (i in 1:nrow(phy$edge)) {
        seg = phy$edge[i, ]
        startDD = meta[seg[1], 2:4]
        stopDD = meta[seg[2], 2:4]
        if (sum(startDD[-3] == stopDD[-3]) == 0) {
            tmp = t(sapply(seq(0, 1, by = 0.05), curvy, startDD, 
                stopDD))
        }
        else {
            tmp = matrix(rep(unlist(startDD[1:2]), each = 3), 
                3, 2)
        }
        pts = cbind(tmp, rep(as.numeric(startDD[3]), nrow(tmp)))
        str = NULL
        for (j in 1:nrow(pts)) {
            str = c(str, paste(pts[j, ], collapse = ",", sep = ""))
        }
        final.coo = meta[seg[2], 2:4]
        final.coo = paste(final.coo, collapse = ",")
        cat("<Placemark>\n<styleUrl>#unselectedLine</styleUrl>\n<LineString>\n<altitudeMode>relativeToGround</altitudeMode>\n<coordinates>\n", 
            paste(str, collapse = "\n"), "\n", final.coo, "\n", 
            "</coordinates>\n</LineString>\n</Placemark>\n", 
            file = goo, append = TRUE)
    }
    cat("</Folder></Document>\n</kml>", file = goo, append = TRUE)
}

bpec.geoTree <- function(bpecout, file="GoogleEarthTree.kml") {
    Output = list()
    writeLines("Creating GoogleEarth Tree plot...")
    treeEdges = bpecout$tree$treeEdges
    clustprob = bpecout$clust$clusterProbs
    count = bpecout$preproc$count
    coordsLocs = bpecout$input$coordsLocs
    dims = dim(bpecout$clust$sampleMeans)[1]
####################################################################
                                        # required libraries igraph, R2G2, ape
####################################################################
                                        # include network.to.newick.r - function
####################################################################
                                        #source("network.to.newick.mod.r")
                                        #source("~/Desktop/bpec-Rsources/network.to.newick_igraph.r")
####################################################################
                                        # load needed library "igraph"
                                        # library("igraph")
                                        # make graph from edgelist
                                        #treeEdgesOut = data.matrix(treeEdges[,1:2])
    treeEdgesOut = data.matrix(bpecout$tree$treeEdges[,1:2])
    
    dimnames(treeEdgesOut) = NULL
    graphEdges = graph.edgelist(treeEdgesOut, directed=TRUE)
                                        # name vertex sequence
                                        #V(graphEdges)$name  = paste("h", sep="", V(graphEdges))
    V(graphEdges)$name  = paste(V(graphEdges))
                                        # remove un-connected vertices
    graphEdgesSub = subgraph.edges(graph=graphEdges, eids=1:length(E(graphEdges)), delete.vertices = TRUE)
                                        #graphEdgesSub
    
                                        # preparation for haplotype-graph plotting
    clustprob[clustprob %in% NaN] = NA
    
                                        # make proportions (rounded integers that sum up to 1000)
                                        #roundint = round(clusterprobs * datsiz[1:count])
    roundInt = 1000 * round(clustprob, 3)
    rowMat = split(roundInt, row(roundInt))
    attributes(rowMat) = NULL
    
                                        # create newick string without lengths
                                        #graphEdges.nwk = network.to.newick(graphEdgesSub)
                                        # or
    graphEdges.nwk = explode(graphEdgesSub)
                                        # string manipulation
    graphEdges.nwk = paste("(",strsplit(graphEdges.nwk,"\\;"),");", sep="")
    
    ## input newick string to create a tree
    
    graphEdgesTree = read.newick(text=graphEdges.nwk)
                                        # remove singletons
    graphEdgesTree = collapse.singles(graphEdgesTree)
    
                                        # add branch length of 1, as R2G2 needs branch lengths
                                        # Note: bpec model does not make inferences about divergence times!
    phy = compute.brlen(graphEdgesTree, 1)
                                        # with the Newick string stored in a file,
    
                                        # load needed package R2G2
                                        #library(R2G2)
                                        # note: the order of lon and lat must be reversed!
  #  print(coordsLocs)
  #  print(dims)
  #  print(ncol(coordsLocs))
    if(ncol(coordsLocs) == dims)
        {
            coordsLocs = cbind(coordsLocs, seq(1, nrow(coordsLocs)))
        }
    coordsLocsSingle = array(0,dim = c(sum(bpecout$preproc$noSamples), 3))
    counter = 1
    for(i in 1:nrow(coordsLocs))
        {
            for(j in (dims+1):ncol(coordsLocs))
                {
                    if(is.na(coordsLocs[i, j]) == FALSE & coordsLocs[i, j] >= 0)
                        {
                            coordsLocsSingle[counter, 1] = coordsLocs[i, 1]
                            coordsLocsSingle[counter, 2] = coordsLocs[i, 2]
                            coordsLocsSingle[counter, 3] = coordsLocs[i, j]
                            counter = counter+1
                        }      
                }
        }
   
    geo.0 =  coordsLocsSingle[, c(1:2, dim(coordsLocsSingle)[2])]
    
    colnames(geo.0) = c("lat", "lon", "taxa")
    taxa=data.frame(taxa= as.numeric(graphEdgesTree[["tip.label"]])); taxa

                                       # merge final dataframe
    coordsLocsSub = merge(taxa, geo.0, by="taxa",all.x = T)
    geo = coordsLocsSub[, c("taxa","lon","lat")]
                                       # converting the tip taxa of the haplotype tree and their corresponding
                                        # geographical coordinates into a KML file to be displayed into Google Earth
                                        #Phylo2GE(geo, phy, 0.05, minAlt = 5000, maxAlt = 100000, goo = "googleearthtree.kml")
    bpec.phylo2GE(geo, phy, resol=1, goo = file)   
    Output$Geo = geo
    Output$Phy = phy

    write.tree(phy, file = "phystring.tre")
    write.csv(geo, file = "geotree.csv")
                                       # set shape of vertex attribute according to haplotypes
                                        # select only sampled haplotypes, the rest should be labeled by haplotype number
  

    return(Output)
}

