bpec.treeEdges = function(MCMCout) {
   
  count = MCMCout$count
  rootProbMean = MCMCout$rootProbs[1:MCMCout$count] + MCMCout$rootProbs[(MCMCout$count+1):(MCMCout$count * 2)]
  root = which.max(rootProbMean)
  levels = MCMCout$levels
  datSiz = MCMCout$noSamples
  clado = MCMCout$clado

  seqLabels = MCMCout$seqsFile[MCMCout$seqLabels]
 
  
 edgeList = NULL
 plotheight = sqrt(2)
 levels = levels + 1

 if (length(seqLabels) < count) {
     seqLabels = c(seqLabels, rep(0, count - length(seqLabels)))
   }
 
 counter = max(seqLabels) + 1
 for(i in 1:count) {
     if (seqLabels[i] == 0) {
         seqLabels[i] = counter
         counter = counter + 1
       }
   }
 
 nodeOrder = array(0, count)
 newNodeOrder = array(0, count)
 maxLevel = max(levels)
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
    # print(descendantcounter)
    # print(prevord)
     descendantcounter = 1
     previousOne = -1
     
     if (prevOrd > 0) {
         for(j in 1:prevOrd) {             
             for(l in 1:count) {
                 if ((clado[(nodeOrder[j] - 1) * count + l] == 1 || clado[(l-1) * count + nodeOrder[j]] == 1) && levels[l] == i) {
                      # print(l)
                                        # edgelist output
                     edgeList = rbind(edgeList, data.frame(vert.from = seqLabels[nodeOrder[j]], vert.to = seqLabels[l], level = i, ssize = datSiz[l]))                   
                     
                     previousOne = l
                     newNodeOrder[descendantCounter] = l
                     descendantCounter = descendantCounter + 1
                     next
                   }
               }
           }
       }
   }
 return(edgeList)
}

