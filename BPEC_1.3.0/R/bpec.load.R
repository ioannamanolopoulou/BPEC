bpec.loadSeq = function(seqsFile) {
  rawSeqs = read.nexus.data(seqsFile)
  return(rawSeqs)
}

bpec.loadCoords = function(coordsFile, header = FALSE) {
    coordsLocs = read.table(coordsFile, header = FALSE, fill = TRUE)
    
    if(header == TRUE) {
        names(coordsLocs) = as.character(unlist(as.matrix(coordsLocs[1, ])))
        coordsLocs = coordsLocs[-1, ]
    }
   coordsLocs = apply(coordsLocs, 2, as.numeric)
 return(coordsLocs)
}

