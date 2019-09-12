bpec.contourPlot <- function(bpecout, GoogleEarth = 0, colorCode = c(7,5,6,3,2,8,4,9), mapType = 'plain', mapCentre = NULL, zoom = 6) {
    writeLines("Creating geographical contour plot...")
    
    if (mapType == 'none') {
        Rmaps = 0
    }
    if (mapType == 'plain') {
        Rmaps = 1
    }
    if (mapType == 'google') {
        Rmaps = 2
    }
    if (mapType == 'osm') {
        Rmaps = 3
    }
    
    rootLocs = numeric(3)
    if(length(bpecout$tree$rootLocProbs) == 1)
        {
            rootLocs[1] = 1
            rootLocs[2] = 1
        }
    if(length(bpecout$tree$rootLocProbs) == 2)
        {
            rootLocs[1] = sapply(sort(bpecout$tree$rootLocProbs, index.return = TRUE), `[`, length(bpecout$tree$rootLocProbs) - 1 + 1)[2]
            rootLocs[2] = sapply(sort(bpecout$tree$rootLocProbs, index.return = TRUE), `[`, length(bpecout$tree$rootLocProbs) - 2 + 1)[2]
            rootLocs[3] = rootLocs[2]
        }
     if(length(bpecout$tree$rootLocProbs) > 2)
        {
            rootLocs[1] = sapply(sort(bpecout$tree$rootLocProbs, index.return = TRUE), `[`, length(bpecout$tree$rootLocProbs) - 1 + 1)[2]
            rootLocs[2] = sapply(sort(bpecout$tree$rootLocProbs, index.return = TRUE), `[`, length(bpecout$tree$rootLocProbs) - 2 + 1)[2]
            rootLocs[3] = sapply(sort(bpecout$tree$rootLocProbs, index.return = TRUE), `[`, length(bpecout$tree$rootLocProbs) - 3 + 1)[2]
        }
    coordsLocs = bpecout$input$coordsLocs[,1:2]
    
    NCP = 15
    NCPR = 2 * NCP + 1
    
    meanSamples = bpecout$clust$sampleMeans
    covSamples = bpecout$clust$sampleCovs
    noClusters=dim(meanSamples)[2]
    fullMeanSamples = meanSamples
    fullCovSamples = covSamples
    
    subSeq=seq(1, length(covSamples[1, 1, 1, ]), length.out = 10)
    meanSamples = meanSamples[, , subSeq]
    covSamples = covSamples[, , , subSeq]
    
    fullClust=numeric(noClusters)
    for(i in 1:noClusters) {
        if (length(which(!is.na(meanSamples[1, i, ]))) > 0) {
            fullClust[i] = 1
        } else {
            fullClust[i] = 0
        }
    }
    
    
    if (Rmaps < 2) {
        xtot=seq(min(min(meanSamples[2, fullClust == 1, ]-sqrt(covSamples[2, 2, fullClust == 1, ]), na.rm = TRUE)),max(max(meanSamples[2, fullClust == 1, ]+sqrt(covSamples[2, 2, fullClust == 1, ]), na.rm = TRUE)), length.out = NCPR)    
        ytot=seq(min(min(meanSamples[1, fullClust == 1, ]-sqrt(covSamples[1, 1, fullClust == 1, ]), na.rm = TRUE)), max(max(meanSamples[1, fullClust == 1, ]+sqrt(covSamples[1, 1, fullClust == 1, ]), na.rm = TRUE)), length.out = NCPR)    
        
        plot(xtot, ytot, col='white', xlab='', ylab='', xlim = range(xtot), ylim = range(ytot), axes=FALSE, asp=1)
                                        #plot(xtot,ytot, col='white',xlab='',ylab='',xlim=range(xtot),ylim=range(ytot),axes=FALSE,asp=1)
    }
    if (Rmaps == 2) {
        colnames(coordsLocs) = c('latitude','longitude')
        coordsLocs = as.data.frame(coordsLocs)
        
        xtot=seq(min(min(meanSamples[2, fullClust == 1 ,] - 3 * sqrt(covSamples[2, 2, fullClust == 1, ]), na.rm = TRUE)), max(max(meanSamples[2, fullClust ==1 , ]+3 * sqrt(covSamples[2, 2, fullClust == 1, ]), na.rm = TRUE)), length.out = NCPR)    
        ytot=seq(min(min(meanSamples[1, fullClust == 1, ] - 3 * sqrt(covSamples[1, 1, fullClust == 1, ]), na.rm = TRUE)), max(max(meanSamples[1, fullClust == 1, ] + 3 * sqrt(covSamples[1, 1, fullClust == 1, ]), na.rm = TRUE)), length.out = NCPR)    
        
        style1 = paste('feature:all|element:labels|visibility:off',
            '&style=feature:road|element:all|visibility:off',
            '&style=feature:administrative|element:geometry|visibility:off',
            '&style=feature:transit.line|element:geometry|visibility:off',
            '&style=feature:water|element:geometry.fill|color:0xffffff',
                                        #    '&style=feature:landscape|element:all|colour:0xffffff',
                                        #     '&style=feature:landscape|element:all|lightness:100',
                                        #      '&style=feature:landscape|element:all|saturation:-100',
                                        # '&style=feature:topography|element:all|colour:0xffffff',
                                        #  '&style=feature:elevation|element:all|colour:0xffffff',
                                        #  '&style=feature:landscape.natural.cover|element:all|colour:0xffffff',
                                        #   '&style=feature:landscape.natural.cover|element:all|lightness:80',
                                        #   '&style=feature:landscape.natural.cover|element:all|saturation:10',
                                        # '&style=feature:landscape.natural.terrain|element:all|colour:0xffffff',
                                        #  '&style=feature:landscape.natural.terrain|element:all|lightness:80',
                                        #  '&style=feature:landscape.natural.terrain|element:all|saturation:10',
            sep="")

        if(!is.null(mapCentre)) {
            totalMap =  ggmap::get_googlemap(mapCentre, zoom = zoom, maptype = "terrain", color = "bw", style = style1, messaging = FALSE)      
            
        } else {  
          #     totalMap =  ggmap::get_map(c(min(xtot), min(ytot), max(xtot), max(ytot)), maptype = "terrain", source = "google", messaging = FALSE)
            totalMap =  ggmap::get_googlemap(c((min(xtot) + max(xtot))/2, (min(ytot) + max(ytot))/2), zoom = zoom, maptype = "terrain", color = "bw", style = style1, messaging = FALSE)
        }
        contourMap = ggmap::ggmap(totalMap) 
        
        
    }
    if (Rmaps == 3) {
        colnames(coordsLocs) = c('latitude','longitude')
        coordsLocs = as.data.frame(coordsLocs)
        xtot=seq(min(min(meanSamples[2, fullClust == 1, ] - 3 * sqrt(covSamples[2, 2, fullClust == 1, ]), na.rm = TRUE)), max(max(meanSamples[2, fullClust == 1, ] + 3 * sqrt(covSamples[2, 2, fullClust == 1, ]), na.rm = TRUE)), length.out = NCPR)    
        ytot=seq(min(min(meanSamples[1, fullClust == 1, ] - 3 * sqrt(covSamples[1, 1, fullClust == 1, ]), na.rm = TRUE)), max(max(meanSamples[1, fullClust == 1, ] + 3 * sqrt(covSamples[1, 1, fullClust == 1, ]), na.rm = TRUE)), length.out = NCPR)         

        if(!is.null(mapCentre)) {
            y1 = max(mapCentre[2] + (max(ytot) - min(ytot))/2, max(ytot))
            x1 = min(mapCentre[1] - (max(xtot) - min(xtot))/2, min(xtot))

            y2 = min(mapCentre[2] - (max(ytot) - min(ytot))/2, min(ytot))
            x2 = max(mapCentre[1] + (max(xtot) - min(xtot))/2, max(xtot))
            
            map = OpenStreetMap::openmap(c(lat = y1, lon = x1), c(lat = y2, lon = x2), type = "osm")
        } else {
            map = OpenStreetMap::openmap(c(lat = max(ytot), lon = min(xtot)), c(lat = min(ytot), lon = max(xtot)), type = "osm")
        }
        mapLatLon = OpenStreetMap::openproj(map)
        contourMap = autoplot.OpenStreetMap(mapLatLon)
        
    }
    nv = length(coordsLocs[, 1])
    newCoordsLocs = array(NA, c(nv, 2))
    for(i in 1:nv) {
        newCoordsLocs[i, 1] = coordsLocs[i, 2] 
        newCoordsLocs[i, 2] = coordsLocs[i, 1]
    }
    
    AT=TRUE
    for(i in 1:noClusters) {     
        for(it in 1:dim(meanSamples)[3]) {
            newMeans = array(NA, c(2, noClusters))
            newCovs = array(NA, c(2, 2, noClusters))
            for(j in 1:noClusters) {
                newMeans[1, j] = meanSamples[2, j, it]
                newMeans[2, j] = meanSamples[1, j, it]
                newCovs[1, 1, j] = covSamples[2, 2, j, it]
                newCovs[1, 2, j] = covSamples[2, 1, j, it]
                newCovs[2, 1, j] = covSamples[1, 2, j, it]
                newCovs[2, 2, j] = covSamples[1, 1, j, it]
            }
           
            if (is.na(newCovs[1,1,i])==0)  {
                x = seq(min(min(meanSamples[2, i, ], na.rm = TRUE)) - 2 * sqrt(max(max(covSamples[2, 2, i, ], na.rm = TRUE))), max(max(meanSamples[2, i, ], na.rm = TRUE)) + 2 * sqrt(max(max(covSamples[2, 2, i, ], na.rm = TRUE))), length.out = NCPR)        
                y = seq(min(min(meanSamples[1, i, ], na.rm = TRUE)) - 2 * sqrt(max(max(covSamples[1, 1, i, ], na.rm = TRUE))), max(max(meanSamples[1, i, ], na.rm = TRUE)) + 2 * sqrt(max(max(covSamples[1, 1, i, ], na.rm = TRUE))), length.out = NCPR)

                xy = expand.grid(x,y)
                z = array(dmvnorm(xy, newMeans[, i], newCovs[, , i]),dim=c(length(x),length(y)))
                
              
                colCont = col2rgb(colorCode[i], alpha = FALSE)
                
                conti = contourLines(x, y, z, levels = 0.5 * max(z))
                if (Rmaps < 2) {
                    polygon(conti[[1]]$x, conti[[1]]$y, col=c(rgb(colCont[1] / 255, colCont[2] / 255, colCont[3] / 255, alpha = 0.1)), border = NA, new = AT, lwd = 2, xlab = "", ylab = "")
                }
                if (Rmaps >= 2) {
                    contour.data = as.data.frame(conti[[1]])
                  
                    contourMap = contourMap + ggplot2::geom_polygon(data = contour.data, ggplot2::aes(x = x, y = y), color = NA, fill = c(rgb(colCont[1] / 255, colCont[2] / 255, colCont[3] / 255)), alpha = 0.1)
                }
                AT=FALSE
            }
        }
                                        #cat ("Press [enter] to continue")
                                        #line <- readline()
    }
    for(i in 1:noClusters) {    
        newMeans[1, i] = mean(meanSamples[2, i,], na.rm = TRUE)
        newMeans[2, i] = mean(meanSamples[1, i,], na.rm = TRUE)
        
        newCovs[1, 1, i] = mean(covSamples[2, 2, i,], na.rm = TRUE)
        newCovs[1, 2, i] = mean(covSamples[2, 1, i,], na.rm = TRUE)
        newCovs[2, 1, i] = mean(covSamples[1, 2, i,], na.rm = TRUE)
        newCovs[2, 2, i] = mean(covSamples[1, 1, i,], na.rm = TRUE)        
    }
    x = array(NA, c(NCPR, noClusters))
    y = array(NA, c(NCPR, noClusters))
    z = array(NA, c(length(x[, 1]), length(y[, 1]), noClusters))
    
    for(i in 1:noClusters) {
        if (fullClust[i] > 0) {
            x[, i] = newMeans[1, i] + (( - NCP:NCP) * 6 * sqrt(newCovs[1, 1, i])) / NCPR
            y[, i] = newMeans[2, i] + (( - NCP:NCP) * 6 * sqrt(newCovs[2, 2, i])) / NCPR

            xy = expand.grid(x[, i],y[, i])
            z[, , i] = array(dmvnorm(xy, newMeans[, i], newCovs[, , i]),dim=c(length(x[,i]),length(y[,i])))              
        }
    }

    AT = TRUE
    contourData = apply(z,c(1,2),sum)
    for(i in 1:noClusters) {
        if (fullClust[i] > 0) {
            if (Rmaps < 2) {
                contour(x[, i], y[, i],z[, , i], col = colorCode[i], levels = 0.5 * max(z[, , i]), xlim = c(min(xtot, na.rm=TRUE), max(xtot, na.rm=TRUE)), ylim =c(min(ytot, na.rm=TRUE), max(ytot, na.rm=TRUE)), add = TRUE, lwd = 1, drawlabels = FALSE, axes = FALSE, xlab = "", ylab = "")
            }
            if (Rmaps >= 2) {               
                contiMean = contourLines(x[,i], y[,i], z[,,i], levels = 0.5 * max(z[, , i]))             
                contourMap = contourMap + ggplot2::geom_path(data = as.data.frame(contiMean[[1]]), ggplot2::aes(x = x, y = y),color = colorCode[i]) 
            }
            AT = FALSE
        }
    }           
    
                                        #  points(newCoordsLocs, pch = 16, cex = 0.3)
    if (Rmaps < 2) {
        points(newCoordsLocs, pch = 17, cex = 0.5, col=1)
        
        points(newCoordsLocs[rootLocs[1], 1], newCoordsLocs[rootLocs[1], 2], pch = 17, cex = 1.1, col = 1)
        points(newCoordsLocs[rootLocs[2], 1], newCoordsLocs[rootLocs[2], 2], pch = 17, cex = 0.9, col = 1)
        points(newCoordsLocs[rootLocs[3], 1], newCoordsLocs[rootLocs[3], 2], pch = 17, cex = 0.7, col = 1)
    }
    if (Rmaps >= 2) {
        contourMap = contourMap + ggplot2::geom_point(data = coordsLocs, ggplot2::aes(x = longitude,y = latitude), size=0.4, color = 'black') 
        
        longitude = NULL
        latitude = NULL
        root1 = data.frame(longitude = newCoordsLocs[rootLocs[1], 1],latitude = newCoordsLocs[rootLocs[1], 2])
        root2 = data.frame(longitude = newCoordsLocs[rootLocs[2], 1],latitude = newCoordsLocs[rootLocs[2], 2])
        root3 = data.frame(longitude = newCoordsLocs[rootLocs[3], 1],latitude = newCoordsLocs[rootLocs[3], 2])
        
        
        contourMap = contourMap + ggplot2::geom_point(data = root1, ggplot2::aes(x = longitude,y = latitude), size = 1.5, color = 'black') 
        contourMap = contourMap + ggplot2::geom_point(data = root2, ggplot2::aes(x = longitude,y = latitude), size = 1.3, color = 'black') 
        contourMap = contourMap + ggplot2::geom_point(data = root3, ggplot2::aes(x = longitude,y = latitude), size = 1.1, color = 'black') 
        
        contourMap = contourMap + labs(x="longitude", y="latitude") 
        plot(contourMap)
    }
    if (Rmaps == 1) {
        world(add = TRUE)
    }
    if (GoogleEarth == 1) {
        
                                        # initialize 
        reslist = list()
        for(i in 1:noClusters) {
            if (i > 1) {
                AT = TRUE
            }
            reslist[[i]] = list(x[, i], y[, i], z[, , i]); 
        }
                                        # end of sections copied from contrgen
        
#### this section produces Google Earth output
##### FUNCTION for kmlPoints output ######
                                        #[R-sig-Geo] Preparing KML files for GoogleMaps/Earth:
                                        #Why no kmlPoint() in maptools package?
                                        #Robert J. Hijmans r.hijmans at gmail.com
        kmlPoints = function(obj, kmlFile, kmlname = "") {
            if (kmlname=="") {kmlname = basename(kmlFile)}
            kmlHeader = c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>","<kml xmlns=\"http://earth.google.com/kml/2.2\">", "<Document>",paste("<name>", kmlname, "</name>", sep = ""))
            kml = ""
            obj[,3] = gsub('&', 'and', obj[,3])
            for (p in 1:length(obj[,1])) {
                                        #		kml = append(kml, "<Folder>")
                kml = append(kml, "<Placemark>")
                kml = append(kml, paste("<name>", obj[p,3], "</name>", sep = ""))
                kml = append(kml, paste("<Point><coordinates>", obj[p,1], ',',obj[p,2], ",0</coordinates></Point>", sep = ""))
                kml = append(kml, "</Placemark>")
                                        #		kml = append(kml, "</Folder>")
            }
            kmlFooter = c("</Document>", "</kml>")
            cat(paste(c(kmlHeader, kml, kmlFooter), sep = "", collapse ="\n"), "\n", file = kmlFile, sep = "")
        }
        
#############################################################
##### MODIFIED FUNCTION for kmlRoots output ######
        root.color = "ff00ff00"
        kmlRoots = function(obj2, kmlFile, kmlname = "") {
            if (kmlname=="") {kmlname = basename(kmlFile)}
            kmlHeader = c("<?xml version=\"1.0\" encoding=\"UTF-8\"?>","<kml xmlns=\"http://earth.google.com/kml/2.2\">", "<Document>",paste("<name>", kmlname, "</name>", sep = ""))
            kml = ""
            obj[,3] = gsub('&', 'and', obj[,3])
            for (p in 1:length(obj2[,1])) {
                                        #		kml = append(kml, "<Folder>")
                kml = append(kml, "<Placemark>")
                kml = append(kml, paste("<name>", obj2[p,3], "</name>", sep = ""))
                kml = append(kml, "<Style>")
                kml = append(kml, "<IconStyle>")
                kml = append(kml, paste("<color>",root.color,"</color>", sep = ""))
                kml = append(kml, paste("<scale>", 1.2, "</scale>", sep = ""))
                kml = append(kml, "<Icon>")
                kml = append(kml, paste("<href>","http://maps.google.com/mapfiles/kml/shapes/arrow.png","</href>", sep = ""))
                kml = append(kml, "</Icon>")
                kml = append(kml, "</IconStyle>")
                kml = append(kml, "</Style>")
                kml = append(kml, paste("<Point><coordinates>", obj2[p,1], ',',obj2[p,2], ",0</coordinates></Point>", sep = ""))
                kml = append(kml, "</Placemark>")
                                        #		kml = append(kml, "</Folder>")
            }
            kmlFooter = c("</Document>", "</kml>")
            cat(paste(c(kmlHeader, kml, kmlFooter), sep = "", 
                      collapse ="\n"), "\n", file = kmlFile, sep = "")
        }
        
        
##### export to GoogleEarth GE, load library with drivers:
                                        # - for Windows: give temporary files the extension .kml, 
                                        #then double-click the *.kml files in the temp directory to start GE
        
#### GENERATION OF PLOTS
        ## loop over lists     
        for(i in 1:noClusters)  {
            if (fullClust[i]>0) {
                x = unlist(reslist[[i]][[1]]);
                                        #x
                y = unlist(reslist[[i]][[2]]);
                                        #y
                z = matrix(unlist(reslist[[i]][[3]]), nrow = length(x), ncol = length(y));
                                        # z
                                        #generate contourLines
                l = contourLines(x, y, z, levels = max(z) * c(0.25, 0.5, 0.75)); 
                ## set geographic projection
                ## here only geographic coordinates are allowed 
                llCRS = CRS("+proj=longlat +datum=WGS84 +towgs84=0,0,0")
################################################################################
                                        # points
                newCoordsLocs.sp = SpatialPoints(newCoordsLocs, llCRS)
                newCoordsLocs.spdf  =  SpatialPointsDataFrame(newCoordsLocs, proj4string=llCRS, data.frame(matrix(seq(1:length(newCoordsLocs[, 1])), ncol = 1))) ##### export sampling points to GoogleEarth
                                        #  str(newCoordsLocs.spdf)
                                        # newCoordsLocs.spdf@data
                obj = cbind(coordinates(newCoordsLocs.sp), newCoordsLocs.spdf@data)
                                        # obj
#########
                                        # contours
                ## create empty list to hold lines
                ll = vector("list", length(l))
                ## loop over list of contours to create Lines object
                for (k in 1:length(l))
                    ll[[k]] = Lines(list(Line(cbind(l[[k]]$x, l[[k]]$y))), as.character(k))
#### make those lines Spatial, and generate a SpatialLinesDataFrame
                ll = SpatialLines(ll, proj4string=llCRS)
#### return the length of the lines slot, how many lines are there?
                out  =  sapply(slot(ll, "lines"),  function(x) {kmlLine(x,name = slot(x, "ID"), col = colorCode[i], lwd=1.5,description = paste("contours_[i]",  slot(x, "ID"))) })
##### generates R temporary files
                                        #tf  =  tempfile()
                tf  =  paste("GoogleEarthContour", as.character(i), ".kml", sep = "")
                kmlFile  =  file(tf,  "w")
                cat(kmlLine(kmlname = paste("Contours", i, sep = "_"),  kmldescription = "<i>BPC</i>")$header, file=kmlFile, sep="\n")
                cat(unlist(out["style",]), file = kmlFile, sep="\n")
                cat(unlist(out["content",]), file = kmlFile, sep="\n")
                cat(kmlLine()$footer, file = kmlFile, sep="\n")
                close(kmlFile)
            }
        }
#############
### plot points
                                        #tf  =  tempfile()
        tf  =  "GoogleEarthPoints.kml"
        kmlFile  =  file(tf,  "w")
        kmlPoints(obj, tf)
        
        llCRS = CRS("+proj=longlat +datum=WGS84 +towgs84=0,0,0")
        newCoordsLocs.sp = SpatialPoints(newCoordsLocs[rootLocs, ], llCRS)
        newCoordsLocs.spdf  =  SpatialPointsDataFrame(newCoordsLocs[rootLocs, ], proj4string = llCRS, data.frame(matrix(seq(1:length(newCoordsLocs[rootLocs, 1])), ncol = 1))) ##### export sampling points to GoogleEarth
                                        #  str(newCoordsLocs.spdf)
                                        # newCoordsLocs.spdf@data
        obj2 = cbind(coordinates(newCoordsLocs.sp), newCoordsLocs.spdf@data)
        close(kmlFile)
        tf  =  "GoogleEarthRoots.kml"
        kmlFile  =  file(tf,  "w")
        kmlRoots(obj2, tf)
        close(kmlFile)
#############
    }
    if(Rmaps >= 2) {
        return(contourMap)
    }
}
