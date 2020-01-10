bpec.treePlot <- function(bpecout,colorCode=c(7,5,6,3,2,8,4,9)) {
  writeLines("Creating clustered tree plot...")
  count = bpecout$preproc$count
  rootProbMean = (bpecout$tree$rootProbs[1, ] + bpecout$tree$rootProbs[2, ]) / 2


  root = which.max(rootProbMean)

  levels = bpecout$tree$levels
  datSiz = bpecout$preproc$noSamples
  clado = bpecout$tree$clado
  clusterProbs = bpecout$clust$clusterProbs
  seqLabels = bpecout$preproc$seqsFile[bpecout$preproc$seqLabels]

  if (length(seqLabels) < count) {
     seqLabels = c(seqLabels, rep(0, count - length(seqLabels)))
   }

  plotHeight = sqrt(2)

  maxMig = dim(clusterProbs)[2]-1
  for(j in 1:count) {
      for(i in 1:(maxMig)) {
          for(k in (i+1):(maxMig+1)) {
              clusterProbs[j, i] = clusterProbs[j, i] + clusterProbs[j, k]
            }
        }
    }

  for(i in 1:count) {
      levels[i] = levels[i] + 1
    }

  nodeCoords = array(0, 2 * count)
  nodeOrder = array(0, count)
  newNodeOrder = array(0, count)

  maxLevel = max(levels)
  plot(0, 0, xlim = c(-0.01,1.01), ylim = c(-sqrt(2) * 1.01, 0.01), type = "n", axes = FALSE, xlab = "", ylab = "")
  nodeCoords[2 * (root-1) + 1] = 1/2
  nodeCoords[2 * (root-1) + 2] = 0
  nodeOrder[1] = root
  newNodeOrder[1] = root
  descendantCounter = 2
  descendants = 0

  for(i in 1:maxLevel) {
      prevDescendants = descendants

      descendants = 1
      nodeOrder = newNodeOrder
      for(j in 1:count) {
          if (levels[j] == i) {
              descendants = descendants + 1
            }
        }
      prevOrd = descendantCounter - 1
      descendantCounter = 1
      previousOne = -1

      if (prevOrd > 0) {
          for(j in 1:prevOrd) {
              for(l in 1:count) {
                  if ((clado[(nodeOrder[j] - 1) * count + l] == 1 || clado[(l-1) * count + nodeOrder[j]] == 1) && levels[l] == i)
                    {
                      nodeCoords[(l-1) * 2 + 1] = descendantCounter / descendants
                      nodeCoords[(l-1) * 2 + 2] = -levels[l] / maxLevel * plotHeight
		      print("now")
		      print(nodeOrder)
		      print(j)
		      print((nodeOrder[j]-1)*2+1)
		      print(nodeCoords)
                      print(previousOne)
                      if ((nodeCoords[(nodeOrder[j] - 1) * 2 + 1] < (descendantCounter + 1) / descendants && (previousOne == -1 || nodeCoords[(nodeOrder[j] - 1) * 2 + 1] > nodeCoords[(previousOne - 1) * 2 + 1])) || descendants == 1) {
                          nodeCoords[(l-1) * 2 + 1] = nodeCoords[(nodeOrder[j] - 1) * 2 + 1]
                        }

                      thickness = bpecout$tree$edgeTotalProb[nodeOrder[j], l] * 1.5
                      lines(c(nodeCoords[(l-1) * 2 + 1], nodeCoords[(nodeOrder[j] - 1) * 2 + 1]), c(nodeCoords[(l-1) * 2 + 2], nodeCoords[(nodeOrder[j] - 1) * 2 + 2]), lwd = thickness)

                      previousOne = l
                      newNodeOrder[descendantCounter] = l
                      descendantCounter = descendantCounter + 1
                      if (i == maxLevel) {
                          if (datSiz[l] > 0) {
                              for(w in 1:(maxMig + 1)) {
                                  bpec.circles(x = nodeCoords[(l-1) * 2 + 1], y = nodeCoords[(l-1) * 2 + 2], radius = 0.02 + 0.0005 * datSiz[l], col = colorCode[w], angle = 360 * (clusterProbs[l, w]))
                                }
                              text(nodeCoords[(l-1) * 2 + 1], nodeCoords[(l-1) * 2 + 2], seqLabels[l], cex = 1)
                            }
                          if (datSiz[l] == 0) {
                              bpec.circles(x = nodeCoords[(l-1) * 2 + 1], y = nodeCoords[(l-1) * 2 + 2], radius = 0.003, col = 1)
                            }
                        }

                      next
                    }
                }
#####
              if (datSiz[nodeOrder[j]] > 0) {
                  for(w in 1:(maxMig + 1)) {
                      bpec.circles(x = nodeCoords[(nodeOrder[j] - 1) * 2 + 1], y = nodeCoords[(nodeOrder[j] - 1) * 2 + 2], radius = 0.02 + 0.0005 * datSiz[nodeOrder[j]], col = colorCode[w], angle = 360 * (clusterProbs[nodeOrder[j], w]))
                    }

                  text(nodeCoords[(nodeOrder[j] - 1) * 2 + 1], nodeCoords[(nodeOrder[j] - 1) * 2 + 2], seqLabels[nodeOrder[j]], cex = 1)
                }
              if (datSiz[nodeOrder[j]] == 0) {
                  bpec.circles(x = nodeCoords[(nodeOrder[j] - 1) * 2 + 1], y = nodeCoords[(nodeOrder[j] - 1) * 2 + 2], radius = 0.003, col = 1)
                }
#####
            }
        }
    }

  output = list()
#  writeLines("Creating GoogleEarth Tree plot...")
  treeEdges = bpecout$tree$treeEdges
  clustProb = bpecout$clust$clusterProbs
  count = bpecout$preproc$count
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
                                        #TreeEdgesOut = data.matrix(TreeEdges[,1:2])
  treeEdgesOut = data.matrix(bpecout$tree$treeEdges[, 1:2])

  dimnames(treeEdgesOut) = NULL
  graphEdges = graph.edgelist(treeEdgesOut, directed = TRUE)
                                        # name vertex sequence
                                        #V(GraphEdges)$name  = paste("h", sep="", V(GraphEdges))
  V(graphEdges)$name = paste(V(graphEdges))
                                        # remove un-connected vertices
  graphEdgesSub = subgraph.edges(graph = graphEdges, eids = 1:length(E(graphEdges)), delete.vertices = TRUE)
                                        #GraphEdgesSub
  roundInt = 1000 * round(clustProb, 3)
  rowMat = split(roundInt, row(roundInt))
  attributes(rowMat) = NULL
  vertShape = ifelse(!is.na(roundInt[c(1:count)]) == T, "pie", "none")

                                        #call the pdf writer
                                        # pdf(file="HaplotypeSubgraphNetwork.pdf", paper="a4r", width = 0, height = 0)
 # par(mai = c(0, 0, 1, 0))
  par(mfrow = c(1, 2))
  set.seed(12345)

  plot.igraph(graphEdgesSub, layout = layout.kamada.kawai,
              vertex.shape = vertShape, vertex.pie = rowMat,
              palette = palette(),
              vertex.pie.color = list(colorCode), vertex.pie.lty = 0,
              vertex.size = 6, vertex.label = V(graphEdgesSub)$name,
              vertex.label.cex = 0.6, vertex.label.font = 2,vertex.label.dist = 0,
              edge.color = "black", edge.arrow.size = 0.05)
  ##run the plot
                                        #  dev.off() #close the device
  output$graphEdgesSub = graphEdgesSub

                                        #TreeEdgesOut = data.matrix(TreeEdges[,1:2])
  treeEdgesOut = data.matrix(bpecout$tree$treeEdges[, 1:2])

  dimnames(treeEdgesOut) = NULL
  graphEdges = graph.edgelist(treeEdgesOut, directed=TRUE)
                                        # name vertex sequence
                                        #V(GraphEdges)$name  = paste("h", sep="", V(GraphEdges))
  V(graphEdges)$name = paste(V(graphEdges))
                                        # remove un-connected vertices
  graphEdgesSub = subgraph.edges(graph = graphEdges, eids = 1:length(E(graphEdges)), delete.vertices = TRUE)
                                        #GraphEdgesSub

                                        # preparation for haplotype-graph plotting
  clustProb[clustProb %in% NaN] = NA

                                        # make proportions (rounded integers that sum up to 1000)
                                        #roundint = round(clusterprobs * datsiz[1:count])
  roundInt = 1000 * round(clustProb, 3)
  rowMat = split(roundInt, row(roundInt))
  attributes(rowMat) = NULL

                                        # create newick string without lengths
                                        #GraphEdges.nwk = network.to.newick(GraphEdgesSub)
                                        # or
  graphEdges.nwk = explode(graphEdgesSub)
                                        # string manipulation
  graphEdges.nwk = paste("(",strsplit(graphEdges.nwk,"\\;"),");", sep="")

  ## input newick string to create a tree

  graphEdgesTree = read.newick(text = graphEdges.nwk)
                                        # remove singletons
  graphEdgesTree = collapse.singles(graphEdgesTree)

                                        #   pdf(file="HaplotreeNoSingles.pdf", paper="a4r", width = 0, height = 0)
  plot(graphEdgesTree, type = "c",
       direction = "l", adj = 1,
       show.node.label = T, label.offset = 0.5, srt = 20,
       font = 3, cex = 0.4,
       root.edge = T)
                                        #   dev.off()
  output$graphEdgesTree = graphEdgesTree
  return(output)
}

