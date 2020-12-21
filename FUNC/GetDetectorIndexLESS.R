GetDetectorIndexLESS <- function(habitat.mx = habitat.mx
                                 ,
                                 detectors.xy = detectors.xy
                                 ,
                                 maxDist = maxDist
                                 ,
                                 ResizeFactor = 1
                                 ,
                                 plot.check = TRUE
){
   
   
   ## ==== 1. CREATE THE XY COORDINATES OF THE HABITAT ====
   habitatID <- habCoordsx <- habCoordsy <- habitat.mx 
   if(ResizeFactor == 1){
      from <- 0.5
      ResizeFactor <- 1 
   }else{  
      from <- (ResizeFactor/2) + 0.5
   }
   ## get dimensions matrix 
   dimCoords <- c(ceiling(dim(habCoordsx)[1]/ResizeFactor), ceiling(dim(habCoordsx)[2]/ResizeFactor))
   CoordsMax <- dimCoords*ResizeFactor
   
   habCoordsx <- matrix(rep(seq(from, CoordsMax[2], by=ResizeFactor ), dimCoords[1]),
                        nrow = dimCoords[1], ncol= dimCoords[2], byrow = T )
   habCoordsy <- matrix(rep(seq(from, CoordsMax[1], by=ResizeFactor ), dimCoords[2]),
                        nrow = dimCoords[1], ncol=dimCoords[2], byrow = F )
   habCoordsy <- as.vector(habCoordsy)
   habCoordsx <- as.vector(habCoordsx)
   
   ## coordinates from habitat 
   habCoordsxy <- cbind(habCoordsx, habCoordsy)

   
   ## ====   2. RESCALE THE HABITAT   ====
   r <- raster(habitat.mx)
   if(ResizeFactor>1){ r <- aggregate(r, fact= ResizeFactor)}
   r[r>0] <-1
   # plot(r)
   habitat.mx1 <- as.matrix(r)
   
   # CREATE A HABITAT ID MATRIX 
   habitatID <- habitat.mx1
   habitatID[] <- as.character(habitat.mx1)
   # ONLY GIVE AN ID TO THE HABITAT CELLS==1 
   habitatID[habitat.mx1 == "1"] <- 1:sum(habitat.mx1=="1")
   m <- sapply(habitatID, FUN=as.numeric)
   habitatID <- matrix(m, nrow=dim(habitatID)[1], ncol=dim(habitatID)[2], byrow =F)
   
   ##REMOVE HABITAT CELL COORDINATES THAT ARE NOT HABITAT  
   habCoordsxy  <- habCoordsxy[as.character(as.vector(habitat.mx1))=="1",]
   
   ## ==== 3. DETERMINE DETECTORS THAT ARE WITHIN A CERTAIN DISTANCE FROM EACH HABITAT CELL  ====
   # Determine detector within radius distance from the center of each habitat cell
   detector.index <- apply(habCoordsxy, 1, function(x){
      D <- sqrt((x[1] - detectors.xy[,1])^2 + (x[2] - detectors.xy[,2])^2) 
      which(D< maxDist)
   })
   
   #make sure it always returns a list. 
   if(class(detector.index)=="matrix"){
      detector.index <- lapply(1:dim(detector.index)[2],function(x) detector.index[,x])
   }
   
   ## ==== 4. STORE DETECTOR INDEX IN A MATRIX ====
   # get number of detectors within the radius for each cell
   nDetectorsLESS <- unlist(lapply(detector.index, function(x) length(x)))
   maxNBDets <- max(nDetectorsLESS)
   # store detector index (colums) for each habitat cell (rows)
   detectorIndex <- matrix(0, nrow=length(detector.index), ncol = maxNBDets)
   for(j in 1:length(detector.index)){
      if(length(detector.index[[j]])!=0){
         detectorIndex[j, 1:nDetectorsLESS[j]] <- detector.index[[j]]
      }
   }
   
   ##PLOT CHECK 
   if(plot.check){
      SXY <- habCoordsxy[sample(1:dim(habCoordsxy)[1],size=1),]# + c(jitter(0,factor = 2), jitter(0,factor = 2) )
      sxyID <- habitatID[trunc(SXY[2]/ResizeFactor)+1, trunc(SXY[1]/ResizeFactor)+1]
      #sxyID <- 80
      
      index <- detectorIndex[sxyID,1:nDetectorsLESS[sxyID]]
      plot(habCoordsxy[,2]~habCoordsxy[,1], pch=16, cex=0.1)
      points(habCoordsxy[sxyID,2]~habCoordsxy[sxyID,1], pch=16, cex=0.4, col="orange")
      
      points(detectors.xy[,2]~detectors.xy[,1], pch=16, cex=0.2, col="red")
      points(detectors.xy[index,2]~detectors.xy[index,1], pch=16, cex=0.4, col="blue")
      points(SXY[2]~SXY[1], bg="red", pch=21, cex=1.2)
      
   }
   
   return <-list( habitatID = habitatID,
                  detectorIndex = detectorIndex,
                  nDetectorsLESS = nDetectorsLESS,
                  maxNBDets = maxNBDets,
                  ResizeFactor = ResizeFactor)
   
   return(return)
}
