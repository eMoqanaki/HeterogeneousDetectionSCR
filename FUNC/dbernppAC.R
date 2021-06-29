#' @title Probability Density Function of the Bernoulli Point Process for Activity Center Placement
#' 
#' @description Probability density function of the Bernoulli point process (i.e. a binomial point process 
#' with one point only) used for activity center placement in SCR studies.  
#' 
#' @param x A vector of coordinates of a single spatial point.
#' @param lowerCoords A matrix of lower coordinates of all habitat windows. One row for one window.
#' Each window should be of size 1x1 (after rescaling if necessary). 
#' @param upperCoords A matrix of upper coordinates of all habitat windows. One row for one window.
#' Both \code{lowerCoords} and \code{upperCoords} are only needed for the \code{rbernppAC} function.
#' @param logIntensities A vector of log habitat selection intensities for all habitat windows. 
#' @param logSumIntensity The log of the total habitat selection intensity over all windows. This can be obtained using
#' \code{logIntensities}. Providing this as an argument helps to reduce some calculations when the function is used repeatedly in a loop.
#' @param habitatGrid A matrix of habitat window indices. When the grid has only one row/column, we need to provide 
#' some (artificial) indices to inflate \code{habitatGrid} so that it is OK to use the function in NIMBLE model code.   
#' @param numGridRows The number of rows of the habitat grid. This may be different from the value obtained 
#' by \code{dim(habitatGrid)[1]}.
#' @param numGridCols The number of columns of the habitat grid. This may be different from the value obtained 
#' by \code{dim(habitatGrid)[2]}.
#' @param log If \code{TRUE} then return the log probability density.
#' 
#' @return The (log) probability density of the observation vector \code{x}.
#' 
#' @author Wei Zhang
#' @export
#' 
dbernppAC <- nimbleFunction(
  run = function(
    x               = double(1),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    logIntensities  = double(1),
    logSumIntensity = double(0),
    habitatGrid     = double(2),
    numGridRows     = integer(0),
    numGridCols     = integer(0),
    log             = integer(0, default = 0)
  ) {
    returnType(double(0))
    ## Check if the point falls within the habitat: 
    ## ***: Getting numGridRows and numGridCols using the following code takes some time
    ## and may cause inefficiency if the function is called repeatedly in a loop.
    ## numGridRows <- dim(habitatGrid)[1]
    ## numGridCols <- dim(habitatGrid)[2]
    ## In addition, when the true habitat grid has one row/column we need to inflate it for use in NIMBLE model code. 
    ## In this case, we have problems.
    ## 
    ## Note that we need to rescale the habitat gird to ensure x and y coordinates start from 0
    ## and each window is of size 1x1. So the following code works correctly. 
    if(min(x) < 0 | x[2] >= numGridRows | x[1] >= numGridCols) {
      if(log) return(-Inf) 
      else return(0.0)
    }
    ## Find which window the point x falls within
    windowInd <- habitatGrid[trunc(x[2])+1, trunc(x[1])+1]
    ## windowInd == 0 means this window is not defined as habitat
    if(windowInd == 0) {
        if(log) return(-Inf)
        else return(0.0)
    }
    ## Log probability density 
    logProb <- logIntensities[windowInd] - logSumIntensity
    
    if(log) return(logProb)
    else return(exp(logProb))
  }
)

#' @title Simulation Function of the Bernoulli Point Process for Activity Center Placement
#' 
#' @description Function to simulate data (one spatial point) from the Bernoulli point process for activity center placement in SCR studies.
#' 
#' @param n The number of samples to generate (this should always be set to 1)
#' @param lowerCoords A matrix of lower coordinates of all habitat windows. One row for one window.
#' Each window should be of size 1x1 (after rescaling if necessary). 
#' @param upperCoords A matrix of upper coordinates of all habitat windows. One row for one window.
#' @param logIntensities A vector of log habitat selection intensities for all habitat windows. 
#' @param logSumIntensity The log of the total habitat selection intensity over all windows. 
#' Providing this as an argument helps to reduce some calculations when the \code{dbernppAC} function is used repeatedly in a loop.
#' Not needed for \code{rbernppAC}.
#' @param habitatGrid A matrix of habitat window indices. Not needed for \code{rbernppAC}.
#' @param numGridRows The number of rows of the habitat grid. Not needed for \code{rbernppAC}.
#' @param numGridCols The number of columns of the habitat grid. Not needed for \code{rbernppAC}.
#' 
#' @return A vector of coordinates of the point generated from the Bernoulli point process.
#' 
#' @author Wei Zhang
#' @export
#' 
rbernppAC <- nimbleFunction(
  run = function(
    n               = integer(0),
    lowerCoords     = double(2),
    upperCoords     = double(2),
    logIntensities  = double(1),
    logSumIntensity = double(0),
    habitatGrid     = double(2),
    numGridRows     = integer(0),
    numGridCols     = integer(0)
  ) {
    returnType(double(1)) 
    if(n <= 0) stop("The number of requested samples must be above zero")
    else if(n > 1) print("rbernppAC only allows n = 1; using n = 1")
    ## Simulate window index
    windowInd <- rcat(1, exp(logIntensities))
    numDims <- 2 ## We consider 2D models for now
    ## A uniform distribution is used within the window
    outCoordinates <- lowerCoords[windowInd,] + 
      runif(numDims, 0.0, 1.0) * (upperCoords[windowInd,] - lowerCoords[windowInd,])
    return(outCoordinates)
  }
)

## Register the distribution
registerDistributions(list(
  dbernppAC = list(
    BUGSdist = "dbernppAC(lowerCoords, upperCoords, logIntensities, logSumIntensity, habitatGrid, numGridRows, numGridCols)",
    types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)", "logIntensities = double(1)", 
              "logSumIntensity = double(0)", "habitatGrid = double(2)", "numGridRows = double(0)", "numGridCols = double(0)"),
    pqAvail = FALSE,
    mixedSizes = TRUE)
))
