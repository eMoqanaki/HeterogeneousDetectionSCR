##########################################################################################
#####        SINGLE-SESSION SCR IGNORING SPATIALLY HETEROGENEOUS DETECTION:          #####
#####                                                                                #####
#####               SIMULATION SCRIPT AND MODEL DEFINITION IN NIMBLE                 #####
#####         CONTINUOUS SPATIAL VARIATION IN BASELINE DETECTION PROBABILITY         #####
#####                           POISSON OBSERVATION PROCESS                          #####
#####                  https://doi.org/10.1007/s10980-021-01283-x                    #####
##########################################################################################

## CLEAN
rm(list = ls())
cat("\014")
gc()

## MEMORY
utils::memory.limit(size = 8e+6) ##--in Mb

## REQUIRED LIBRARIES
library(rgeos)              
library(rgdal)              
library(raster)            
library(coda)               
library(nimble)             
library(ggplot2)            
library(R.utils)            
library(abind)              
library(nimbleSCR)

## SOURCE THE REQUIRED FUNCTIONS
dir.function <- paste(getwd(), "/FUNC", sep = "") ##--source FUNC directory on your machine
R.utils::sourceDirectory(dir.function, modifiedOnly = FALSE)

inv.logit <- function(inValues) {
  1.0 / (1.0 + exp(-inValues))
}


## ----------------------------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## ----------------------------------------------------------------------------------------------

### ==== 0.1. GENERAL VARIABLES DECLARATION ====
myVars <- list( 
  
  ## SIMULATION PARAMETERS
  HABITAT =    list(extent     = 20,      # habitat extent
                    resolution = 1,       # habitat resolution
                    buffer     = 4.5),    # buffer around the habitat: 3*sigma
  DETECTORS =  list(resolution = 1),      # detector resolution 
  DENSITY =    list(Beta.dens  = 0),      # set to 0 for uniform location of ACs
  POPULATION = list(N = 250,              # true population size
                    M = 625),             # size of the augmented population: 2.5*N
  DETECTIONS = list(p0 = 0.15,            # baseline detection probability (intercept)   
                    det.betas.p0 = -0.5,  # a factor covariate with 2 levels
                    sigma = 1.5,          # the scale parameter
                    phi   = 1000,         # defines the level of spatial autocorrelation: High (0.001), Intermediate (1), Low (1000)
                    det.prop1 = 1)        # the proportion of detectors with higher detectability
)


## ----------------------------------------------------------------------------------------------
## ------ 1. SET-UP HABITAT AND DETECTORS -----
## ----------------------------------------------------------------------------------------------

### ==== 1.1. GENERATE HABITAT ====
buffer <- myVars$HABITAT$buffer
grid.size <- c(myVars$HABITAT$extent)

coords <- matrix(c(buffer            , buffer ,
                   grid.size + buffer, buffer ,
                   grid.size + buffer, grid.size + buffer,
                   buffer            , grid.size + buffer,
                   buffer            , buffer 
), ncol = 2, byrow = TRUE)

P1 <- Polygon(coords)
myStudyArea <- SpatialPolygons(list(Polygons(list(P1), ID = "a")),
                               proj4string = CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

## MAKE A BUFFER
buff <- rgeos::gBuffer(myStudyArea, width = myVars$HABITAT$buffer)
myStudyArea.buff <- raster::union(myStudyArea, buff)
r <- raster(extent(myStudyArea.buff)) ##--Rasterise polygon
res(r) <- myVars$HABITAT$resolution   ##--Change resolution
r[] <- rep(0,ncell(r))

habitat.r <- rasterize(as(myStudyArea.buff, "SpatialLines"), r, field=1) ##--rasterise as lines and then polygons so all cells touched by the polygon are covered
habitat.r <- rasterize(myStudyArea.buff, habitat.r, update=T, field=1)
habitat.r[is.na(habitat.r)] <- 0
habitat.mx <- as.matrix(habitat.r)

IDCells.r <- habitat.r
IDCells.r[] <- 1:length(IDCells.r) ##--Cell ID starts from the top left corner
IDCells.mx <- as.matrix(IDCells.r)

## Obtain xy coordinates of cells
habitat.xy <- xyFromCell(habitat.r, 1:ncell(habitat.r))
dimnames(habitat.xy) <- list(1:length(habitat.xy[,"x"]), c("x","y"))
habitat.sp <- SpatialPointsDataFrame(data.frame(habitat.xy[,c("x","y")]), data = data.frame(habitat.xy), proj4string = CRS(projection(myStudyArea)))

## Obtain lower and upper cell coordinates
lowerHabCoords <- data.frame(coordinates(habitat.sp) - myVars$HABITAT$resolution/2)
upperHabCoords <- data.frame(coordinates(habitat.sp) + myVars$HABITAT$resolution/2)

nHabCells <- dim(lowerHabCoords)[1]
habQualityDims <- dim(habitat.r)[2:1]

### ==== 1.2. GENERATE DETECTORS ====
det.r <- rasterize(myStudyArea, habitat.r, field = 1) ##--DEFINE A MATRIX WITHOUT THE BUFFER AREA TO CALCULATE N EXCLUDING BUFFER WITHIN THE NIMBLE MODEL

detector.r1 <- raster(extent(myStudyArea), resolution = 1, crs = proj4string(myStudyArea))
detector.r1[] <- rep(1,length(detector.r1))

main.detector.xy <- xyFromCell(detector.r1, 1:ncell(detector.r1))
colnames(main.detector.xy) <- c("x","y")
main.detector.sp <- SpatialPointsDataFrame(data.frame(main.detector.xy[,c("x","y")]),
                                           data=data.frame(main.detector.xy),
                                           proj4string=CRS(projection(myStudyArea)))

plot(habitat.r)
plot(myStudyArea, add = T)
points(main.detector.sp)

### ==== 1.3. SCALE DETECTORS & HABITAT & UPPER/LOWER COORDINATES FOR HABITAT WINDOWS ====
scaled <- scaleCoordsToHabitatGrid(coordsData = main.detector.xy,
                                   coordsHabitatGridCenter = habitat.xy)

ScaledLowerCoords <- scaleCoordsToHabitatGrid(coordsData = lowerHabCoords,
                                              coordsHabitatGridCenter = habitat.xy)$coordsDataScaled
ScaledUpperCoords <- scaleCoordsToHabitatGrid(coordsData = upperHabCoords,
                                              coordsHabitatGridCenter = habitat.xy)$coordsDataScaled
ScaledUpperCoords[,2] <- ScaledUpperCoords[,2]+1
ScaledLowerCoords[,2] <- ScaledLowerCoords[,2]-1

## CREATE A HABITAT GRID
habIDCells.mx <- IDCells.mx
habIDCells.mx[] <- 0
scaledHabGridCenters <- scaled$coordsHabitatGridCenterScaled[habitat.r[]==1,]

for(i in 1:nrow(scaledHabGridCenters)){
  habIDCells.mx[trunc(scaledHabGridCenters[i,2])+1,
                trunc(scaledHabGridCenters[i,1])+1] <- i
}#i

### ==== 1.4. CREATE CACHED OBJECTS FOR LESS RESTRICTION ====
DetectorIndexLESS <- getLocalObjects(habitatMask = habitat.mx,
                                     coords = scaled$coordsDataScaled,
                                     dmax = 10/res(habitat.r)[1],##--this value should be large enough so that all detections fall withing the restricted window
                                     resizeFactor = 1,
                                     plot.check = F)


## ----------------------------------------------------------------------------------------------
## ------ 2. DEFINE THE PARAMETER SPACE -----
## ----------------------------------------------------------------------------------------------

#### ==== 2.1. DETECTOR COVARIATE ====
## simulate autocorrelated rasters
MyDensityRaster <- SimulateAutocorrelatedRaster(sp = main.detector.sp,
                                                NRaster = 1,
                                                phi = myVars$DETECTIONS$phi,
                                                scaled = TRUE,
                                                focal = NULL,
                                                plot = F)


par(mfrow = c(2,2))
plot(MyDensityRaster)
hist(MyDensityRaster[])

## ******************************************************************************************************************
## TURN COVARIATE INTO A UNIFORMLY DISTRIBUTED ONE (instead of some random normal, with SE varying according to Moran's index)
df <- data.frame(row = 1:length(MyDensityRaster[]), orig.value = (MyDensityRaster[])[,1])
df <- df[order(df$orig.value),]

df$value <- seq(-1.96, 1.96, length.out = length(df$orig.value)) ##--UNIFORM
df <- df[order(df$row), ]

MyDensityRaster[] <- df$value
plot(MyDensityRaster)
hist(MyDensityRaster[])
## ******************************************************************************************************************

## extract detector covs for a continuous cov
main.detector.sp$detCov <- data.frame(raster::extract(MyDensityRaster, main.detector.sp))[1]
colnames(main.detector.sp$detCov) <- "detCov"

### ==== 2.2. SIMULATE INDIVIDUAL AC LOCATIONS ====
## Sample random (uniform) activity center locations
simulated.ACS <- sp::spsample(x = raster::aggregate(raster::rasterToPolygons(habitat.r,function(x)x>0)),
                              n = myVars$POPULATION$N,
                              type = "random")

simulated.ACS$id <- 1:length(simulated.ACS)

### ==== 2.3. SIMULATE DETECTION ARRAY : y.ar ====
## EXPORT EXAMPLE STATE SPACE AND DETECTION FIGURES
y <- array(0, c(length(simulated.ACS), length(main.detector.sp)))
D <- rgeos::gDistance(main.detector.sp, simulated.ACS, byid = TRUE)

temp <- formula(paste("~", paste(names(main.detector.sp$detCov), collapse = "+"), "-1", sep = ""))
Xmat <- model.matrix(temp, main.detector.sp$detCov)
det.fixed.effects.p0 <- Xmat%*%myVars$DETECTIONS$det.betas.p0
det.fixed.effects.p0 <- do.call(rbind, lapply(1:length(simulated.ACS), function(x)t(det.fixed.effects.p0)))
X.det.p0 <- data.frame(Xmat)

intercept.p0 <- matrix(log(myVars$DETECTIONS$p0), length(simulated.ACS), length(main.detector.sp))
p0 <- exp(intercept.p0 + det.fixed.effects.p0)
p <- p0 * exp(-D * D/(2 * myVars$DETECTIONS$sigma * myVars$DETECTIONS$sigma))

y[] <- apply(p, c(1, 2), function(x) rpois(1, x))

### ==== 2.4. SUBSET TO INDIVIDUALS DETECTED AT LEAST ONE YEAR/DATA SET ==== 
## CHECK THE NUMBER OF INIDVIDUALS DETECTED
detected <- apply(y,1, function(x) sum(x)>0)
sum(detected)

## SUBSET UNDETECTED IDS
y <- y[detected,]

### ==== 2.5. AUGMENT DATA ====

## Augmentation factor that ensures that the total number of individuals 
## (alive + available) is always the same, regardless of the simulation
this.AF <- myVars$POPULATION$M/sum(detected)-1

sum(detected) * (1+ this.AF) == myVars$POPULATION$M

y <- rbind(y, matrix(0, nrow = myVars$POPULATION$M-sum(detected), ncol = dim(y)[2]))
z <- c(rep(1, sum(detected)), rep(NA,myVars$POPULATION$M-sum(detected)))

### ==== 2.6. TRANSFORM Y TO SPARSE MATRICES =====

SparseY <- getSparseY(y)

### ==== 2.7.   SET THE INPUT FOR NIMBLE ====
### ==== 2.7.1. Define the nimble model  ====

modelCode <- nimbleCode({
  
  sumIntensity <- sum(mu[1:numHabWindows])
  
  ##---- SPATIAL PROCESS 
  for(i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = mu[1:numHabWindows],
      logSumIntensity = sumIntensity,
      numGridRows = y.max,
      numGridCols = x.max,
      habitatGrid = habitatGrid[1:y.max,1:x.max])
  }#i
  
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  
  for (i in 1:n.individuals){
    z[i] ~ dbern(psi)
  }#i
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0, 50)
  alpha <- -1 / (2 * sigma * sigma)
  p0 ~ dunif(0, 10)
  
  for (i in 1:n.individuals){
    
    ## LESS APPROACH AND USE OF SPARSE MATRICES
    y[i, 1:nMaxDetectors] ~ dpoisLocal_normal(detNums = nbDetections[i],
                                              detIndices = yDets[i,1:nMaxDetectors],
                                              lambda = p0,
                                              s = sxy[i,1:2],
                                              sigma = sigma,
                                              trapCoords = detector.xy[1:n.detectors,1:2],
                                              localTrapsIndices = detectorIndex[1:n.cells,1:maxNBDets],
                                              localTrapsNum = nDetectorsLESS[1:n.cells],
                                              resizeFactor = ResizeFactor,
                                              habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                              indicator = z[i])
  }#i
  
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
})


### ==== 2.7.2. Define the constants to be used in the model ====
## Set the list of model constants
nimConstants <- list(
  x.max = habQualityDims[1],
  y.max =  habQualityDims[2],
  n.detectors = dim(y)[2], 
  n.individuals = dim(y)[1],
  numHabWindows = max(habIDCells.mx),
  nMaxDetectors = SparseY$maxDetNums,
  y.maxDet = dim(DetectorIndexLESS$habitatGrid)[1],
  x.maxDet = dim(DetectorIndexLESS$habitatGrid)[2],
  ResizeFactor = DetectorIndexLESS$resizeFactor,
  maxNBDets = DetectorIndexLESS$numLocalIndicesMax,
  n.cells = dim(DetectorIndexLESS$localIndices)[1]
)


### ==== 2.7.3. Define the data to be used in the model ====

## Set the list of data
nimData <- list(
  lowerHabCoords = ScaledLowerCoords,
  upperHabCoords = ScaledUpperCoords,
  z = z,
  y = SparseY$y[,,1],
  detector.xy = scaled$coordsDataScaled,
  mu = rep(1, max(habIDCells.mx)),
  yDets = SparseY$detIndices[,,1],
  nbDetections = SparseY$detNums[,1],
  detectorIndex = DetectorIndexLESS$localIndices,
  nDetectorsLESS = DetectorIndexLESS$numLocalIndices,
  habitatIDDet = DetectorIndexLESS$habitatGrid,
  habitatGrid = habIDCells.mx
)

### ==== 2.7.4. Define the initial values to be used in the model  ====
## Initialise AC locations
sxy <- MakeInitXY(y = y,
                  habitat.mx = habitat.mx,
                  detector.xy = scaled$coordsDataScaled,
                  IDCells.mx = IDCells.mx,
                  grid.xy = scaled$coordsHabitatGridCenterScaled)

## Initialise z values
z.init <- ifelse(!is.na(z), NA, 1)
z.init[!is.na(z.init)] <- rbinom(sum(!is.na(z.init)), size = 1, prob = 0.5)


nimInits <- list(sigma = runif(1,0,10),
                 p0 = runif(1,0,10),   
                 z = z.init,           
                 psi = runif(1,0,1),
                 sxy = sxy[ , ,1]      
)

### ==== 2.7.5. Define the parameters to be monitored ====
nimParams <- c("N","sxy","z","sigma","p0","psi")


## ----------------------------------------------------------------------------------------------
## ------ 3. SCR - A TEST RUN
## ----------------------------------------------------------------------------------------------

### ==== 3.1. CREATE & COMPILE THE NIMBLE MODEL OBJECT & MCMC CONFIGURATION ====
model <- nimbleModel(code = modelCode,
                     constants = nimConstants,
                     data = nimData,
                     inits = nimInits,
                     check = FALSE,       
                     calculate = FALSE)  

cmodel <- compileNimble(model)
cmodel$calculate()

MCMCconf <- configureMCMC(model = model,
                          monitors = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1) 

MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)

### ==== 3.2. RUN THE MCMC ====
myNimbleOutput <- runMCMC(mcmc = cMCMC,
                          nburnin = 1, 
                          niter = 1000,#15000  
                          nchains = 3, 
                          samplesAsCodaMCMC = TRUE)


coda::traceplot(myNimbleOutput[,"N"])
coda::traceplot(myNimbleOutput[,"sigma"])
coda::traceplot(myNimbleOutput[,"p0"])