#' @title Function to calculate new coordinates from a SpatialPointsDataFrame object based on a another SpatialPointsDataFrame that contains locations of grid cell centers of a regular grid.
#' #'
#' @description
#' \code{UTMToGrid} returns a dataframe with x and y coordinates.
#' 
#' @param data.sp \code{SpatialPointsDataframe} with with points for which new coordinates are to be calculated. 
#' @param grid.sp \code{SpatialPointsDataframe} with with points at grid cell centers  of the grid that will form the basis of the transformation to from utm to grid units.
#' @param plot.check Logical for whether \code{True} or not (\code{FALSE}) plots are to be generated execution.
#' @param data.sxy An \code{array} up to 5 dimensions with the posterior sxy coordinates. xy are assumed to be placed on the 3 dimension of the array.  
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/UTMToGrid.R
#' @keywords SCR data prep
#'
#' @examples
#' # UTMToGrid:
#' 
#'UTMToGrid(data.sp = detector.sp, grid.sp = habitat.sp, plot.check = TRUE)
#'UTMToGrid(data.sp = habitat.sp, grid.sp = habitat.sp, plot.check = TRUE)
#
#' 
#' 
#'  

UTMToGrid <- function(
   data.sp,
   grid.sp,
   data.sxy = NULL,
   plot.check = TRUE)
{
   # PREPARE THE DATA
   grid.xy <- as.array(coordinates(grid.sp))
   dimnames(grid.xy) <- list(1:length(grid.sp), c("x","y"))
   
   # CHECK IF THE DATA ARE THE POSTERIOR SXY 
   if(is.null(data.sxy)){
      data.xy <- as.array(coordinates(data.sp))
      dimnames(data.xy) <- list(1:length(data.sp), c("x","y"))
   }else{
      data.xy <- data.sxy
   }
   
   # CALCULATE THE RESOLUTION
   resolution <- min(diff(unique(sort(grid.xy[ ,"x"]))))#---assumes square grid cells and utm projection (units in meters or km; not latlong!)
   
   
   ## obtain x and y min
   start0.y <- max(grid.xy[ ,"y"]) + resolution/2 #---because we are moving from top to bottom
   start0.x <- min(grid.xy[ ,"x"]) - resolution/2 #---because we are moving from left to right
   
   ##---- TO CHECK: re-projecting the grid cell centers
   grid.scaled.xy <- grid.xy
   
   grid.scaled.xy[ ,"y"] <- (start0.y - grid.xy[ ,"y"])/resolution
   grid.scaled.xy[ ,"x"] <- (grid.xy[ ,"x"] - start0.x)/resolution
   
   ##---- REPROJECTING THE DATA
   if(length(dim(data.xy))==2){
      data.scaled.xy <- data.xy
      data.scaled.xy[ ,"y"] <- (start0.y - data.xy[ ,"y"])/resolution
      data.scaled.xy[ ,"x"] <- (data.xy[ ,"x"] - start0.x)/resolution 
   }
   ##  DEAL WITH CASES WHERE WE HAVE SXY 
   # ITEATIONS,i,xy
   if(length(dim(data.xy))==3){
      data.scaled.xy <- data.xy
      data.scaled.xy[ ,2,] <- (start0.y - data.xy[ ,2,])/resolution
      data.scaled.xy[ ,1,] <- (data.xy[ ,1,] - start0.x)/resolution
   }
   # ITEATIONS,i,xy,t
   if(length(dim(data.xy))==4){
      data.scaled.xy <- data.xy
      data.scaled.xy[ ,,2,] <- (start0.y - data.xy[,,2,])/resolution 
      data.scaled.xy[ ,,1,] <- (data.xy[ ,,1,] - start0.x)/resolution 
   }
   # ITEATIONS,i,xy,t,??
   if(length(dim(data.xy))==5){
      data.scaled.xy <- data.xy
      data.scaled.xy[ ,2,,,] <- (start0.y - data.xy[,2,,,])/resolution 
      data.scaled.xy[ ,1,,,] <- (data.xy[ ,1,,,] - start0.x)/resolution 
   }
   
   ##---- TO CHECK: re-projecting the grid cell centers
   
   if(is.null(data.sxy)){
      if(plot.check==TRUE){
         if(all(par()$mfrow==c(1,1))) par(mfrow=c(1,2), mar=c(5,5,1,1))
         
         plot(rbind(grid.sp), col="white", main="Original")
         plot(grid.sp, add=TRUE, pch=19, col="darkgreen", cex=0.6)
         plot(data.sp, pch=19, col="orange", add=TRUE)
         plot((0-y) ~ x,rbind(data.scaled.xy, grid.scaled.xy), type="n", main="Gridded", axes=FALSE)
         
         axis(1)
         axis(2, at=seq(0, 0 - max(round(data.scaled.xy[ ,"x"],0)), -ceiling(0.1*max(round(data.scaled.xy[ ,"x"], 0)))), 
              labels=seq(0, max(round(data.scaled.xy[,"x"],0)), ceiling(0.1 * max(round(data.scaled.xy[ ,"x"], 0)))))
         points((0-y) ~ x, grid.scaled.xy, pch=19, col="darkgreen", cex=0.6)
         temp.y <- 0-(1:max(round(grid.scaled.xy[ ,"y"], 0)))
         temp.x <- 1:max(round(grid.scaled.xy[ ,"x"],0))
         segments(-10000, temp.y, 10000, temp.y, col="darkgreen")
         segments(temp.x, -10000, temp.x, 10000, col="darkgreen")
         points((0-y)~x, data.scaled.xy, pch=19, col="orange")
      }
   }else{print("Plotting function not available for data.sxy")}
   
   out <- list(grid.scaled.xy = grid.scaled.xy,
               grid.xy = grid.xy,
               data.scaled.xy = data.scaled.xy,
               data.xy = data.xy)
   
   return(out)
}
