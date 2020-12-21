GetSparseY <- function( y = y
                       , Nodetections = -1
){
   # IF SNAPSHOT, CONVERT TO ARRAY
   if(length(dim(y))==2){
      Y <- array(y, c(dim(y),1))
   }else{Y <- y}
   
   # NUMBER OF DETECTIONS FOR EACH ID
   nbDetections <- apply(Y, c(1,3), function(x) length(which(x>0)))

   ySparseDets <- array(-1, c(dim(Y)[1],max(nbDetections), dim(Y)[3]))
   ySparse <- array(-1, c(dim(Y)[1],max(nbDetections), dim(Y)[3]))
   
   # FILL IN THE ARRAYS
   for(t in 1:dim(Y)[3]){
    for(i in 1:dim(Y)[1]){
      if(nbDetections[i,t]>0){
         # GET WHERE (DETECTOR ID) DETECTIONS OCCUR
         ySparseDets[i, 1:nbDetections[i,t],t] <- which(Y[i,,t]>0)
         # GET NUMBE OF DETECTIONS 
         ySparse[i, 1:nbDetections[i,t],t] <- Y[i, which(Y[i,,t]>0),t]
      }
    }
   }   
  

   return(list( y = ySparse ## Detection array 
               ,yDets = ySparseDets
               ,nbDetections = nbDetections
               ,nMaxDetectors = max(nbDetections)))
   
}