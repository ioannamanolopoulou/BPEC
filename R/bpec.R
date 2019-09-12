`bpec` <- function(seqsFile, coordsFile, dims = 2, iter = 100000, postSamples = 100, maxMig = 5, ds = 0, colorCode = c(7,5,6,3,2,8)) {

  rawSeqs = read.nexus.data(seqsFile)
  coordsLocs = bpec.loadCoords(coordsFile)
  
  #run the Markov chain Monte Carlo sampler
  bpecout = bpec.mcmc(rawSeqs, coordsLocs, maxMig, iter, ds, postSamples, dims)
  return(bpecout)
}

`plot.bpec` <-  function(x, GoogleEarth = 0, colorCode = c(7,5,6,3,2,8), ... ) {
  # plot geographical cluster contour map
    bpec.contourPlot(x, GoogleEarth = 0, colorCode = colorCode) 
    
                                        # plot tree network with cluster indicators
    bpec.Tree = bpec.treePlot(x, colorCode = colorCode)
    
                                        # now also plot the environmental covariates
    par(mfrow = c(2, 3)) #split the plot window into 2x3 to fit all the covariates
    if (x$input$coordsDims > 2) {
        bpec.covariatesPlot(x, colorCode = colorCode) 
    }
    if (GoogleEarth == 1) {
        bpec.Geo = bpec.geoTree(x,file="GoogleEarthTree.kml")
    }
}

`summary.bpec` <- function(object,...) {
  cat(paste('bpec ran with the following settings:\n'))
  cat(paste('Number of input sequences',object$input$seqCountOrig,'\n'))
  cat(paste('Input sequence length',object$input$seqLengthOrig,'\n'))
  cat(paste('Number of iterations',object$input$iter,'\n'))
  cat(paste('Number of saved iterations',dim(object$clust$sampleMeans)[3]/2,'\n'))
  
  cat(paste('Dimensions',object$input$coordsDims,'\n'))
  cat(paste('Parsimony relaxation',object$input$ds,'\n'))
  cat(paste('Maximum number of migrations',length(object$tree$migProbs) - 1,'\n'))
  
  cat(paste('\nThe results of bpec are: \n'))
  cat(paste('Number of haplotypes (including missing)',dim(object$preproc$seq)[1],'\n'))
  cat(paste('of which',dim(object$preproc$seq)[1]-object$input$seqCountOrig,'are missing\n'))
  cat(paste('Effective sequence length',object$preproc$seqLength,'\n'))
  cat(paste('The most likely number of migrations is',which.max(object$tree$migProbs)))
  seqLabels = object$preproc$seqsFile[object$preproc$seqLabels]
  maxMig = length(object$tree$migProbs) - 1
  rootProbMean = (object$tree$rootProbs[1, ] + object$tree$rootProbs[2, ]) / 2
 
  if (which.max(rootProbMean) <= length(seqLabels)) {
    root = seqLabels[which.max(rootProbMean)]
    writeLines(paste("\nThe most likely root node is ", root, sep=""))
  }
  if (which.max(rootProbMean) > length(seqLabels)) {
    root = which.max(rootProbMean)
    writeLines(paste("\nThe most likely root node is extinct", sep=""))
  }
  maxVar = 0
  chainMeans = array(0, dim = c(dim(object$clust$sampleMeans)[1], dim(object$clust$sampleMeans)[2], 2))
  flag = 0
  postSamples = length(object$clust$sampleMeans[1, 1, ]) / 2
  
  for(i in 1:(maxMig+1)) {
    for(j in 1:object$input$coordsDims) {
      for(l in 1:2) {
        chainMeans[j, i, l] = mean(object$clust$sampleMeans[j, i, (1 + (l-1) * postSamples):(l * postSamples)], na.rm = TRUE)
      }
      if (sum(is.na(chainMeans[j, i, ])) == 0) {
        if (abs(chainMeans[j, i, 1] - chainMeans[j, i, 2]) > 0.05 * (max(object$input$coordsLocs[, j]) - min(object$input$coordsLocs[, j])) && mean(is.na(object$clust$sampleMeans[j, i, ])) < 0.5) {
        writeLines("NO CLUSTER CONVERGENCE: You need to re-run the sampler with more iterations")
          flag = 1
          break
        }
      }
    }
    if (flag == 1) {
      break
    }        
  }
  
  rootProbs1 = object$tree$rootProbs[1:object$preproc$count]
  rootProbs2 = object$tree$rootProbs[(object$preproc$count + 1):(object$preproc$count * 2)]
  
  if (max(abs(rootProbs1 / sum(rootProbs1) - rootProbs2 / sum(rootProbs2))) > 0.5 / object$preproc$count) {
    writeLines("NO ROOT CONVERGENCE: You need to re-run the sampler with more iterations");
  }   
}


`print.bpec` <- function(x,...) {  
  cat(paste('\nThe results of bpec are: \n'))
  cat(paste('Number of haplotypes (including missing)',dim(x$preproc$seq)[1],'\n'))
  cat(paste('of which',dim(x$preproc$seq)[1]-x$input$seqCountOrig,'are missing\n'))
  cat(paste('Effective sequence length',x$preproc$seqLength,'\n'))
  cat(paste('The most likely number of migrations is',which.max(x$tree$migProbs)))
  seqLabels = x$preproc$seqsFile[x$preproc$seqLabels]
  maxMig = length(x$tree$migProbs) - 1
  rootProbMean = (x$tree$rootProbs[1, ] + x$tree$rootProbs[2, ]) / 2

  if (which.max(rootProbMean) <= length(seqLabels)) {
    root = seqLabels[which.max(rootProbMean)]
    writeLines(paste("\nThe most likely root node is ", root, sep=""))
  }
  if (which.max(rootProbMean) > length(seqLabels)) {
    root = which.max(rootProbMean)
    writeLines(paste("\nThe most likely root node is extinct", sep=""))
  }
  maxVar = 0
  chainMeans = array(0, dim = c(dim(x$clust$sampleMeans)[1], dim(x$clust$sampleMeans)[2], 2))
  flag = 0
  postSamples = length(x$clust$sampleMeans[1, 1, ]) / 2
  
  for(i in 1:(maxMig+1)) {
    for(j in 1:x$input$coordsDims) {
      for(l in 1:2) {
        chainMeans[j, i, l] = mean(x$clust$sampleMeans[j, i, (1 + (l-1) * postSamples):(l * postSamples)], na.rm = TRUE)
      }
      if (sum(is.na(chainMeans[j, i, ])) == 0) {
       if (abs(chainMeans[j, i, 1] - chainMeans[j, i, 2]) > 0.05 * (max(x$input$coordsLocs[, j]) - min(x$input$coordsLocs[, j])) && mean(is.na(x$clust$sampleMeans[j, i, ])) < 0.5) {
         writeLines("NO CLUSTER CONVERGENCE: You need to re-run the sampler with more iterations")
          flag = 1
          break
        }
      }
    }
    if (flag == 1) {
      break
    }        
  }
  
  rootProbs1 = x$tree$rootProbs[1:x$preproc$count]
  rootProbs2 = x$tree$rootProbs[(x$preproc$count + 1):(x$preproc$count * 2)]
  
  if (max(abs(rootProbs1 / sum(rootProbs1) - rootProbs2 / sum(rootProbs2))) > 0.5 / x$preproc$count) {
    writeLines("NO ROOT CONVERGENCE: You need to re-run the sampler with more iterations");
  }   
}

`mean.bpec` <- function(x, ...) {
    cat('The mean cluster centres are:\n')
    print(apply(x$clust$sampleMeans, c(1,2), 'mean', na.rm = TRUE))

    cat('The mean cluster covariances are:\n')
    print(apply(x$clust$sampleCovs, c(1,2,3), 'mean', na.rm = TRUE))

    cat('The mean root posterior probabilities are:\n')
    print(apply(x$tree$rootProbs, 2, 'mean',  na.rm = TRUE))
}

input <- function(bpecout) UseMethod("input")

`input.bpec` <- function(bpecout) {
    output = list()
    output$seqCountOrig = bpecout$input$seqCountOrig
    output$seqLengthOrig = bpecout$input$seqLengthOrig
    output$iter = bpecout$input$iter
    output$ds = bpecout$input$ds
    output$coordsLocs = bpecout$input$coordsLocs
    output$coordsDims = bpecout$input$coordsDims
    output$locNo = bpecout$input$locNo
    output$locData = bpecout$input$locData
 
 return(output)
}

preproc <- function(bpecout) UseMethod("preproc")
`preproc.bpec` <- function(bpecout){
    output = list()
    output$seq = bpecout$preproc$seq
    output$seqsFile = bpecout$preproc$seqsFile
    output$seqLabels = bpecout$preproc$seqLabels
    output$seqIndices = bpecout$preproc$seqIndices
    output$seqLength = bpecout$preproc$seqLength
    output$noSamples = bpecout$preproc$noSamples
    output$count = bpecout$preproc$count

    return(output)
}

#setGeneric("tree", function(bpecout) standardGeneric("tree"))
output.tree <- function(bpecout) UseMethod("output.tree")

`output.tree.bpec` <- function(bpecout){
    output = list()
    output$clado = bpecout$tree$clado
    output$levels = bpecout$tree$levels
    output$edgeTotalProb = bpecout$tree$edgeTotalProb
    output$rootProbs = bpecout$tree$rootProbs
    output$treeEdges = bpecout$tree$treeEdges
    output$rootLocProbs = bpecout$tree$rootLocProbs
    output$migProbs = bpecout$tree$migProbs
    
    return(output)
}

output.clust <- function(bpecout) UseMethod("output.clust")

`output.clust.bpec` <- function(bpecout){
    output = list()
    output$sampleMeans = bpecout$clust$sampleMeans
    output$sampleCovs = bpecout$clust$sampleCovs
    output$sampleIndices = bpecout$clust$sampleIndices
    output$clusterProbs = bpecout$clust$clusterProbs
  
    return(output)
}

#setGeneric("mcmc", function(x) standardGeneric("mcmc"))
#mcmc <- function(x) UseMethod("mcmc")
#setMethod("mcmc", "bpec",
#function(x, ...) {
#  output = list()
#    output$MCMCparams = bpecout$mcmc$MCMCparams
#    output$MCChainMeans = bpecout$mcmc$MCChainMeans
#   output$codaInput = bpecout$mcmc$codaInput 
#    
#    return(output)
#  }
#)
#
#setMethod("mcmc", "list",
#function(x, ...) {
#    mcmc.list(x)
#  }
#)

#setGeneric("mcmc", function(x) standardGeneric("mcmc"))

output.mcmc <- function(bpecout) UseMethod("output.mcmc")

`output.mcmc.bpec` <- function(bpecout){
    output = list()
    output$MCMCparams = bpecout$mcmc$MCMCparams
    output$MCChainMeans = bpecout$mcmc$MCChainMeans
    output$codaInput = bpecout$mcmc$codaInput 
    
    return(output)
}

#`mcmc.list` <- function(bpecout){
#    return(coda::mcmc(bpecout))
#}
