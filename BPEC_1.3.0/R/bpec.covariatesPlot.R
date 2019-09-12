bpec.covariatesPlot = function(bpecout, colorCode = c(7,5,6,3,2,8,4,9)) {
  writeLines("Creating covariate distribution plot...")
  covNames = colnames(bpecout$input$coordsLocs) 
  meanSamples = bpecout$clust$sampleMeans
  covSamples = bpecout$clust$sampleCovs
  noClusters=dim(meanSamples)[2]
  dims=dim(meanSamples)[1]
  subSeq=seq(1, length(covSamples[1, 1, 1, ]), length.out = 20)
  meanSamplesSub = meanSamples[, , subSeq]
  
  fullClust = numeric(noClusters)
  for(i in 1:noClusters) {
    if (length(which(!is.na(meanSamplesSub[1, i, ]))) > 0) {
      fullClust[i] = 1
    } else {
      fullClust[i] = 0
    }
  }
 
  if (dims > 2) {      
    for(i in 3:dims) {
      # w=seq(min(Means[i,])-2*max(sqrt(Covs[i,i,])),max(Means[i,])+2*max(sqrt(Covs[i,i,])),length=100)
      # w=seq(min(MeanSamples[i,,],na.rm=TRUE)-2*max(sqrt(CovSamples[i,i,,]),na.rm=TRUE),max(MeanSamples[i,,],na.rm=TRUE)+2*max(sqrt(CovSamples[i,i,,]),na.rm=TRUE),length=100)
      # w=seq(min(MeanSamples[i,,],na.rm=TRUE)-1*mean(sqrt(CovSamples[i,i,,]),na.rm=TRUE),max(MeanSamples[i,,],na.rm=TRUE)+1*mean(sqrt(CovSamples[i,i,,]),na.rm=TRUE),length=100)
      w=seq(mean(meanSamples[i, , ], na.rm = TRUE) - 3 * mean(sqrt(covSamples[i, i, , ]), na.rm = TRUE), mean(meanSamples[i, , ], na.rm = TRUE) + 3 * mean(sqrt(covSamples[i, i, , ]), na.rm = TRUE), length = 100)
      maxax = 0
      
      for(j in 1:noClusters) {
        if (fullClust[j] == 0) {
          next
        }                
        for(it in 1:dim(meanSamples)[3]) {  
          if (is.na(meanSamples[i, j, it]) == FALSE) {
            dnst = dnorm(w, meanSamples[i, j, it], sqrt(covSamples[i, i, j, it]))
            if (max(dnst) > maxax) {
              maxax = max(dnst)              
            }                 
          }
        }           
      }
      
      plot(w, dnst, col='white', ylim = c(0, maxax), type = 'l', xlab = covNames[i], ylab = "")     
      
      for(j in 1:noClusters) {
        if (fullClust[j] == 0) {
          next
        }
        densIt = array(NA, dim = c(length(w), dim(meanSamples)[3]))
        for(it in 1:dim(meanSamples)[3])  {                      
          if (is.na(meanSamples[i, j, it]) == FALSE) {
            densIt[, it] = dnorm(w, meanSamples[i, j, it], sqrt(covSamples[i, i, j, it]))                        
          }
        }                              
        densUp = w
        densLow = w
        densMed = w
        for(jj in 1:length(w)) {
          densUp[jj] = quantile(densIt[jj, ], 0.9, na.rm = TRUE)
          densLow[jj] = quantile(densIt[jj, ], 0.1, na.rm = TRUE)
          densMed[jj] = quantile(densIt[jj, ], 0.5, na.rm = TRUE)
        }
        colCont = col2rgb(colorCode[j], alpha = FALSE)
        polygon(c(w, rev(w)), c(densLow, rev(densUp)), col = c(rgb(colCont[1] / 255, colCont[2] / 255, colCont[3] / 255, alpha = 0.3)), border = NA)
        lines(w, densMed, col = colorCode[j])        
      }        
    }
  }
}
