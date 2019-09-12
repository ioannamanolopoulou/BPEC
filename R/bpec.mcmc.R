bpec.mcmc <- function(rawSeqs, coordsLocs, maxMig, iter, ds, postSamples=0, dims=-1)
{
   #   cord.dec = SpatialPoints(cbind(coordsLocs[,2], -coordsLocs[,1]), proj4string = CRS("+proj=longlat"))
   #   cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32748"))
   #   coordsLocs[,1:2] = cord.UTM@coords

   #   cord.longlat <- spTransform(cord.UTM, CRS("+proj=longlat +datum=WGS84"))
   #   cord.UTM <- cord.UTM@coords
   #   cord.longlat <- cord.longlat@coords
    
  if (dims == -1){
    dims = round(length(which((round(coordsLocs, digits = 0) == coordsLocs) == FALSE)) / nrow(coordsLocs), digits = 0)
  }
  if (dims !=- 1) {
    dimscal=round(length(which((round(coordsLocs, digits = 0) == coordsLocs) == FALSE)) / nrow(coordsLocs), digits = 0)
    if (dims != dimscal) {
      writeLines("\nThe dataset appears to have some integer-valued environmental/phenotypic entries.\nIf this is not correct, check the value you have given for dims.")
    }
  }
  if (postSamples > iter / 10) {
    writeLines("You cannot have postSamples greater than iter/10, so it will be changed to iter/10.")
    postSamples = iter / 10
  }
  if (iter %% 10 != 0) {
    writeLines("iter needs to be a multiple of 10. The MCMC will not run. ")
    return(-1)
  }
  seqNames = as.numeric(gsub("[^0-9]", "", names(rawSeqs))  )
  seqCount = as.numeric(length(rawSeqs))  
  seqLength = as.numeric(length(rawSeqs[[1]]))
  rawSeqsOrig = rawSeqs
  
  coordinates = coordsLocs[ , 1:dims]
  coordsLocs[is.na(coordsLocs)] = -1
  if (ncol(coordsLocs) == dims) {
    coordsLocs = cbind(coordsLocs, seq(1, nrow(coordsLocs)))
    writeLines("The dataset did not contain a list of haplotypes per row,\nso the program will assume that row 1 corresponds to haplotype 1, row 2 to haplotype 2 etc. ")
  }
  endRow = 1:nrow(coordsLocs)
  endRow[] = -2
  coords = coordsLocs[, 1:2]
  coordsLocs = cbind(coordsLocs, endRow)
  locNo = as.numeric(nrow(coordsLocs))
  maxLoc = as.numeric(ncol(coordsLocs))

  if(sum(coordsLocs[, dims + 1] %in% seqNames) < length(coordsLocs[, dims + 1] )) {
      missingSeq = which(coordsLocs[, dims + 1] %in% seqNames == FALSE)
     #   print(coordsLocs[, dims + 1] %in% seqNames)
     #   print(missingSeq)
      writeLines("Some of the sequences in the location file does not appear in the NEXUS file. The MCMC will not run.")
      return(-1)
  }


  coordsLocsOrig = coordsLocs
  coordsLocs = as.numeric(as.vector(t(coordsLocs)))
  
  #  UnSeqs=RawSeqs
  #  UnCount=SeqCount

   
  rawSeqsList = unlist(rawSeqs)
  rawSeqsList = matrix(rawSeqsList, nrow = length(rawSeqs[[1]]), ncol = length(rawSeqs))
  rawSeqsList = t(rawSeqsList)
  
  delCol = NULL
  for(i in 1:ncol(rawSeqsList)) {
    if (length(unique(rawSeqsList[, i])) == 1) {    
      delCol = c(delCol, i) 
    }  
  }
  rm(rawSeqsList)
  if (is.null(delCol) == FALSE) {
    for(i in 1:length(rawSeqs)) {
      rawSeqs[[i]] = rawSeqs[[i]][-delCol]
    }
  }
  
 
  seqLength = as.numeric(length(rawSeqs[[1]]))

  if(seqLength == 0) {
      for(i in 1:seqCount) {          
          rawSeqs[[i]][1] = 'A'
      }
      seqLength = 1
  }
  seqs = array(NA, dim = (seqCount * seqLength))
  counter = 1
  for (i in 1:seqCount) {
    for(j in 1:seqLength) {
            seqs[counter] = rawSeqs[[i]][j]
            counter = counter+1
        }
  }
  
  count = seqCount
  seeds = 2
  ancestral = 1

  MCMCout = .C("bpecfunction", modeInitial = 1, seqR = seqs, coordsLocsR = coordsLocs, coordsDimsR = as.numeric(dims), seqsFileR = seqNames, seqCountR = count, seqLengthR = seqLength, locNoR = locNo, maxLocR = maxLoc, maxMigR = maxMig, seedsR = seeds, iterR = iter, dsR = ds, ancestralR = ancestral, locDataR = numeric((length(coordsLocs) - dims * locNo) * dims), sampleMeansR = numeric(2), sampleCovsR = numeric(2), sampleIndicesR = numeric(2), sampleClusterCodaR = numeric(2), sampleRootCodaR = numeric(2), postSamplesR = postSamples, levelsR = numeric(2), cladoR = numeric(2), edgeTotalProbR = numeric(2), noSamplesR = numeric(10 * count + 100), clusterProbsR = numeric(2), countR = numeric(2), migPMigProbsR = numeric(maxMig + 1), rootProbsR = numeric(2), rootLocProbsR = numeric(locNo), MCMCparamsR = numeric(8), seqLabelsR = numeric(10 * count + 100), nSeqR = numeric(10), errorCodeR = numeric(1), PACKAGE = "BPEC")
 
  MCMCout$countR = max(c(count, MCMCout$countR[1]))

  if (MCMCout$errorCodeR[1]!=1) {
    MCMCout = .C("bpecfunction", modeInitial = 0, seqR = seqs, coordsLocsR = coordsLocs, coordsDimsR = as.numeric(dims), seqsFileR = seqNames, seqCountR = seqCount, seqLengthR = seqLength, locNoR = locNo, maxLocR = maxLoc, maxMigR = maxMig, seedsR = seeds, iterR = iter, dsR = ds, ancestralR = ancestral, locDataR = numeric((length(coordsLocs) - dims * locNo) * dims), sampleMeansR = numeric((postSamples+1) * dims * (maxMig+1) * seeds), sampleCovsR = numeric((postSamples) * dims * dims * (maxMig+1) * seeds), sampleIndicesR = numeric(sum((MCMCout$noSamplesR+1)) * 2 * (postSamples)), sampleClusterCodaR = numeric((postSamples) * (dims + 4 + (dims-2)) * (maxMig + 1) * seeds), sampleRootCodaR = numeric(postSamples * seeds), postSamplesR = postSamples, levelsR = numeric(MCMCout$countR+1), cladoR = numeric((MCMCout$countR + 1) * (MCMCout$countR + 1)), edgeTotalProbR = numeric((MCMCout$countR + 1) * (MCMCout$countR + 1)), noSamplesR = MCMCout$noSamplesR, clusterProbsR = numeric((MCMCout$countR + 1) * (maxMig + 1)), countR = MCMCout$countR, migProbsR = numeric(maxMig+1), rootProbsR = numeric(2 * MCMCout$countR), rootLocProbsR = numeric(locNo), MCMCparamsR = numeric(8), seqLabelsR = numeric(MCMCout$countR + 1), nSeqR = as.numeric(MCMCout$countR), errorCodeR = numeric(1), PACKAGE = "BPEC")

    names(MCMCout$MCMCparamsR) = c('PsiPrior', 'CentralSieve', 'TrialAdd', 'AveClusterWeight', 'ClustAccRate', 'RootAccRate', 'MeanGamma', 'MeanW')
    
    MCMCout$coordsLocsR = t(array(MCMCout$coordsLocsR, dim = c(length(MCMCout$coordsLocsR) / MCMCout$locNoR, MCMCout$locNoR)))
    MCMCout$coordsLocsR = MCMCout$coordsLocsR[, -ncol(MCMCout$coordsLocsR)]
    
    MCMCout$clusterProbsR = t(array(MCMCout$clusterProbsR, dim = c(maxMig+1, MCMCout$countR+1)))
    MCMCout$countR = MCMCout$countR[1]
    
    #  MCMCout$rootfreqsR = MCMCout$rootfreqsR[1:MCMCout$countR]
    #MCMCout$clustersR=MCMCout$clustersR[1:sum(MCMCout$NoSamplesR)]
    MCMCout$levelsR = MCMCout$levelsR[1:MCMCout$countR]
    MCMCout$clusterProbsR = MCMCout$clusterProbsR[1:MCMCout$countR, ]
    MCMCout$cladoR = MCMCout$cladoR[1:(MCMCout$countR * MCMCout$countR)]
    MCMCout$edgeTotalProbR = MCMCout$edgeTotalProbR[1:(MCMCout$countR * MCMCout$countR)]
    MCMCout$edgeTotalProbR = array(MCMCout$edgeTotalProbR, dim = c(MCMCout$countR, MCMCout$countR))
    MCMCout$locDataR = matrix(MCMCout$locDataR, ncol = dims, byrow = TRUE )
    
    #  MCMCout$RootProbsR = MCMCout$RootProbsR[1:MCMCout$countR]
    
    if (length(MCMCout$seqsFile) < max(c(MCMCout$seqLabelsR, MCMCout$seqLabelsR))) {
      MCMCout$seqsFile=c(MCMCout$seqsFile, rep(0, max(c(MCMCout$seqLabelsR, MCMCout$seqLabelsR))-length(MCMCout$seqsFile)))
    }
    
    MCMCout$sampleMeansR = is.finite(MCMCout$sampleMeansR) * MCMCout$sampleMeansR
    MCMCout$sampleCovsR = is.finite(MCMCout$sampleCovsR) * MCMCout$sampleCovsR
  
    MCMCout$sampleMeansR = array(MCMCout$sampleMeansR, dim=c(dims, maxMig+1, seeds * postSamples))
    MCMCout$sampleClusterCodaR = array(MCMCout$sampleClusterCodaR, dim=c(dims+2 + dims, maxMig+1, seeds * postSamples))
    MCMCout$sampleIndicesR = array(MCMCout$sampleIndicesR, dim = c(sum(MCMCout$noSamplesR), seeds * postSamples))
    
    MCMCout$sampleCovsR = array(MCMCout$sampleCovsR, dim=c(dims, dims, maxMig+1, seeds * postSamples))
    
    maxVar = 0
    chainMeans = array(0, dim = c(dim(MCMCout$sampleMeansR)[1], dim(MCMCout$sampleMeansR)[2],2))
    clusterNonConvergence = numeric(maxMig+1)
    flag = 0
    for(i in 1:(maxMig + 1)) {
      for(j in 1:dims) {
        for(l in 1:2) {
          chainMeans[j,i,l] = mean(MCMCout$sampleMeansR[j, i, (1+(l-1) * postSamples) : (l * postSamples)], na.rm=TRUE)
        }
        if (sum(is.na(chainMeans[j, i, ])) == 0) {    
          if (abs(chainMeans[j, i, 1] - chainMeans[j, i, 2]) > 0.05 * (max(coordinates[, j]) - min(coordinates[, j])) && mean(is.na(MCMCout$sampleMeansR[j, i, ])) < 0.5) {
              clusterNonConvergence[i] = 1
              if (flag == 0) {
                      writeLines("NO CLUSTER CONVERGENCE: You need to re-run the sampler with more iterations")
                      flag = 1
                  }
          }
        }
      }         
    }

    

    rootProbs1 = MCMCout$rootProbsR[1:MCMCout$countR]
    rootProbs2 = MCMCout$rootProbsR[(MCMCout$countR+1):(MCMCout$countR * 2)]
    # print(rootprobs1/sum(rootprobs1))
    # print(rootprobs2/sum(rootprobs2))
    
    if (max(abs(rootProbs1 / sum(rootProbs1)-rootProbs2 / sum(rootProbs2)))>0.5 / MCMCout$countR) {
      writeLines("NO ROOT CONVERGENCE: You need to re-run the sampler with more iterations");
    }
    
    MCMCout$rootProbMeanR = MCMCout$rootProbsR[1:MCMCout$countR] + MCMCout$rootProbsR[(MCMCout$countR+1):(MCMCout$countR * 2)]
    MCMCout$treeEdgesR = bpec.treeEdges(MCMCout)
    
    seqLabels = MCMCout$seqsFileR[MCMCout$seqLabelsR]
    if (which.max(MCMCout$rootProbMeanR) <= length(seqLabels)) {
      root = seqLabels[which.max(MCMCout$rootProbMeanR)]
      writeLines(paste("The most likely root node is ", root, sep=""))
    }
    if (which.max(MCMCout$rootProbMeanR) > length(seqLabels)) {
      root = which.max(MCMCout$rootProbMeanR)
      writeLines(paste("The most likely root node is extinct", sep=""))
    }
    MCMCout$rootProbMeanR = MCMCout$rootProbMeanR / sum(MCMCout$rootProbMeanR)
    
    uniqueLocs = coordinates[!duplicated(coordinates[1:2]), ]
    rootProbsUnique = MCMCout$rootLocProbsR
    MCMCout$locNoR = length(uniqueLocs)
    
    for (i in 1:length(uniqueLocs[, 1])) {          
      identicalLocs = (coordinates[, 1] == uniqueLocs[i, 1]) & (coordinates[, 2] == uniqueLocs[i, 2])
      tempRoot = sum(rootProbsUnique[identicalLocs])
      rootProbsUnique[identicalLocs] = tempRoot
      if (sum(identicalLocs) > 1) {
        firstLoc = min(which(identicalLocs == TRUE))
        identicalLocs[firstLoc] = FALSE
        rootProbsUnique[identicalLocs] = 0
      }      
    }
    
    MCMCout$rootLocProbsR = rootProbsUnique
    
    rootlocs = numeric(3)
    rootlocs[1] = sapply(sort(MCMCout$rootLocProbs, index.return=TRUE), `[`, length(MCMCout$rootLocProbs)-1+1)[2]
    rootlocs[2] = sapply(sort(MCMCout$rootLocProbs, index.return=TRUE), `[`, length(MCMCout$rootLocProbs)-2+1)[2]
    rootlocs[3] = sapply(sort(MCMCout$rootLocProbs, index.return=TRUE), `[`, length(MCMCout$rootLocProbs)-3+1)[2]
    writeLines(paste("The most likely ancestral locations are ",rootlocs[1], ",",rootlocs[2], ",",rootlocs[3], sep=""))
    
    haploIndex = numeric(sum(MCMCout$noSamplesR))
    counter = 1
    for(i in 1:length(MCMCout$noSamplesR)) {
      if (MCMCout$noSamplesR[i] > 0) {
        haploIndex[counter:(counter+MCMCout$noSamplesR[i] - 1)] = i
        counter = counter + MCMCout$noSamplesR[i]
      }
    }

    codatemp = MCMCout$clusterProbsR
    if(MCMCout$countR>1) {        
        for(i in 1:MCMCout$countR) {
            for(j in 1:MCMCout$maxMigR) {
                codatemp[i, j] = sum(MCMCout$sampleIndicesR[haploIndex == i , ] == j)
            }      
            codatemp[i, ] = codatemp[i, ] / sum(codatemp[i, ])
        }
    } else {
        for(j in 1:MCMCout$maxMigR) {
            codatemp[j] = sum(MCMCout$sampleIndicesR[haploIndex == 1 , ] == j)
        }      
        codatemp = codatemp / sum(codatemp)
    }
    
    MCMCout$seqR = MCMCout$seqR[1:(MCMCout$seqLengthR * length(rawSeqs))]
    MCMCout$seqR = array(MCMCout$seqR, dim = c(length(rawSeqs), MCMCout$seqLengthR))
    
    MCMCout$clusterProbsR = codatemp
    # print(MCMCout$ClusterProbsR)
    
    #  if (coda == 1)
    #  {
    dim1 = length(MCMCout$sampleClusterCodaR[, 1, 1])
    dim2 = length(MCMCout$sampleClusterCodaR[1, , 1])
    dim3 = length(MCMCout$sampleClusterCodaR[1, 1, ])
    chain1 = array(MCMCout$sampleClusterCodaR[,,1:(dim3 / 2)], dim=c(dim1*dim2, dim3/2))
   
    chain1 = rbind(chain1,MCMCout$sampleRootCodaR[1 : (dim3 / 2)])
    chain2 = array(MCMCout$sampleClusterCodaR[, , (dim3 / 2 + 1) : (dim3)], dim = c(dim1 * dim2, dim3 / 2))
    chain2 = rbind(chain2,MCMCout$sampleRootCodaR[(dim3 / 2 + 1) : (dim3)])
    MCMCout$codaInput = list()
    MCMCout$codaInput$line1 = coda::mcmc(t(chain1), start = 1, end = dim3/2, thin = 1)
    MCMCout$codaInput$line2 = coda::mcmc(t(chain2), start = 1, end = dim3/2, thin = 1)
    MCMCout$codaInput = coda::mcmc.list(MCMCout$codaInput)
    MCMCout$seqLengthOrig = as.numeric(length(rawSeqsOrig[[1]]))
    MCMCout$seqCountOrig = as.numeric(length(rawSeqsOrig))
    # }
    
    MCMCout$sampleClusterCodaR = NULL
    MCMCout$modeInitial = NULL
    MCMCout$ancestralR = NULL
    MCMCout$maxMigR = NULL
    MCMCout$codaR = NULL
    MCMCout$postSamplesR = NULL
  #   MCMCout$iterR = NULL
   #  MCMCout$dsR = NULL
    MCMCout$nSeqR = NULL
    MCMCout$maxLocR = NULL
    MCMCout$seedsR = NULL
    MCMCout$rootR = NULL
    MCMCout$rootlocsR = NULL
    MCMCout$clustersR = NULL
    MCMCout$seqCorrR = NULL
    
  }

 
  bpecout = list()
  
  bpecout$input = list()
  bpecout$input$seqCountOrig = MCMCout$seqCountOrig
  bpecout$input$seqLengthOrig = MCMCout$seqLengthOrig
  bpecout$input$iter = MCMCout$iterR
  bpecout$input$ds = MCMCout$dsR
  bpecout$input$coordsLocs = MCMCout$coordsLocsR
  colnames(bpecout$input$coordsLocs) = 1:dim(bpecout$input$coordsLocs)[2]
  colnames(bpecout$input$coordsLocs)[1:MCMCout$coordsDimsR] = colnames(coordsLocsOrig)[1:MCMCout$coordsDimsR]
 
  bpecout$input$coordsDims = MCMCout$coordsDimsR
  bpecout$input$locNo = MCMCout$locNoR
  bpecout$input$locData = MCMCout$locDataR
  bpecout$input$locData = bpecout$input$locData[as.logical(rowSums(bpecout$input$locData != 0)), ]
  
  bpecout$preproc = list()
  bpecout$preproc$seq = MCMCout$seqR
 
  bpecout$preproc$seqsFile = MCMCout$seqsFileR
  bpecout$preproc$seqLabels = MCMCout$seqLabelsR
  bpecout$preproc$seqIndices = lapply(rawSeqs, function(elemnt) which(as.logical(lapply(unique(rawSeqs),function(uniquelement) all(uniquelement == elemnt)))))
  bpecout$preproc$seqLength = MCMCout$seqLengthR
  bpecout$preproc$noSamples = MCMCout$noSamplesR
  bpecout$preproc$count = MCMCout$countR
  
  
  bpecout$clust = list()
  bpecout$clust$sampleMeans = MCMCout$sampleMeansR
 
  rownames(bpecout$clust$sampleMeans) = colnames(coordsLocsOrig)[1:dim(bpecout$clust$sampleMeans)[1]]
  rownames(bpecout$clust$sampleMeans)[1:2] = c('lat','lon')
  bpecout$clust$sampleCovs = MCMCout$sampleCovsR
 
  dimnames(bpecout$clust$sampleCovs)[[1]] = colnames(coordsLocsOrig)[1:dim(bpecout$clust$sampleCovs)[1]]
  dimnames(bpecout$clust$sampleCovs)[[1]][1:2] = c('lat','lon')

  dimnames(bpecout$clust$sampleCovs)[[2]] = colnames(coordsLocsOrig)[1:dim(bpecout$clust$sampleCovs)[2]]
  dimnames(bpecout$clust$sampleCovs)[[2]][1:2] = c('lat','lon')
  bpecout$clust$sampleIndices = MCMCout$sampleIndicesR
  bpecout$clust$clusterProbs = MCMCout$clusterProbsR
  
  bpecout$tree = list()
  bpecout$tree$clado = MCMCout$cladoR
  templength = length(bpecout$preproc$seqsFile[bpecout$preproc$seqLabels])
   bpecout$tree$levels = MCMCout$levelsR
  names(bpecout$tree$levels)[1:templength] = bpecout$preproc$seqsFile[bpecout$preproc$seqLabels]
  bpecout$tree$edgeTotalProb = MCMCout$edgeTotalProbR
  rownames(bpecout$tree$edgeTotalProb) = 1:dim(bpecout$tree$edgeTotalProb)[1]
  colnames(bpecout$tree$edgeTotalProb) = 1:dim(bpecout$tree$edgeTotalProb)[2]
  rownames(bpecout$tree$edgeTotalProb)[1:templength] = bpecout$preproc$seqsFile[bpecout$preproc$seqLabels]
  colnames(bpecout$tree$edgeTotalProb)[1:templength] = bpecout$preproc$seqsFile[bpecout$preproc$seqLabels]
  bpecout$tree$rootProbs = t(array(MCMCout$rootProbsR[1:(2*bpecout$preproc$count)], dim = c(bpecout$preproc$count,2)))
  
  bpecout$tree$rootProbs[1, ] = bpecout$tree$rootProbs[1, ] / sum(bpecout$tree$rootProbs[1, ])
  bpecout$tree$rootProbs[2, ] = bpecout$tree$rootProbs[2, ] / sum(bpecout$tree$rootProbs[2, ])
  
  colnames(bpecout$tree$rootProbs) = 1:dim(bpecout$tree$rootProbs)[2]
  colnames(bpecout$tree$rootProbs)[1:templength] = bpecout$preproc$seqsFile[bpecout$preproc$seqLabels]
 
  bpecout$tree$treeEdges = MCMCout$treeEdgesR
  
  bpecout$tree$rootLocProbs = MCMCout$rootLocProbsR
  bpecout$tree$migProbs = MCMCout$migProbsR

  bpecout$mcmc = list()
 
  bpecout$mcmc$MCMCparams =MCMCout$MCMCparamsR
  bpecout$mcmc$codaInput = MCMCout$codaInput
 
  bpecout$mcmc$errorCode = MCMCout$errorCodeR
     
  attr(bpecout, "class") <- "bpec"
  return(bpecout)

 # attr(MCMCout, "class") <- "bpec"
 # return(MCMCout)
}

