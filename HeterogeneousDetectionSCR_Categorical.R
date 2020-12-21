################################################################################
#####        SINGLE-SESSION SCR IGNORING HETEROGENEOUS DETECTION:          #####
#####          SIMULATION SCRIPT AND MODEL DEFINITION IN NIMBLE            #####
#####                                                                      #####
#####       CONTINUOUS VARIATION IN BASELINE DETECTION PROBABILITY         #####
########################################################################### #####


##-- REQUIRED LIBRARIES
library(rgeos)              
library(rgdal)              
library(raster)            
library(coda)               
library(nimble)             
library(ggplot2)            
library(R.utils)            
library(abind)              


##-- SOURCE THE REQUIRED FUNCTIONS
##-- SOURCE THE REQUIRED FUNCTIONS
dir.function <- ## SOURCE FUNC DIRECTORY ON YOUR MACHINE
sourceDirectory(dir.function, modifiedOnly = FALSE)



## ----------------------------------------------------------------------------------------------
## ------ 0. SET ANALYSIS CHARACTERISTICS -----
## ----------------------------------------------------------------------------------------------

### ==== 0.1. GENERAL VARIABLES DECLARATION ====
myVars <- list( 
  
  ## SIMULATION PARAMETERS
  HABITAT =    list(extent     = 20,            # habitat extent
                    resolution = 1,             # habitat resolution
                    buffer     = 4.5),          # buffer around the habitat: 3*sigma
  DETECTORS =  list(resolution = 1),            # detector resolution 
  DENSITY =    list(Beta.dens  = 0),            # set to 0 for uniform location of ACs
  POPULATION = list(N = 250,                    # true population size
                    M = 625),                   # size of the augmented population: 2.5*N
  DETECTIONS = list(p0 = 0.15,                  # baseline detection probability   
                    det.betas.p0 = c(-10000),   # a factor covariate with 2 levels
                    sigma = 1.5,                # the scale parameter
                    phi   = 1000,               # defines the level of spatial autocorrelation
                    det.prop1 = 0.25)           # the proportion of detectors with higher detectability
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

myHabitat <- MakeHabitat(poly = myStudyArea,
                         resolution = myVars$HABITAT$resolution, 
                         buffer = myVars$HABITAT$buffer,
                         plot.check = FALSE)          


lowerHabCoords <- coordinates(myHabitat$lower.hab.sp)
upperHabCoords <- coordinates(myHabitat$upper.hab.sp)
nHabCells <- dim(lowerHabCoords)[1]
habQualityDims <- dim(myHabitat$habitat.r)[2:1]


## DEFINE A MATRIX WITHOUT THE BUFFER AREA TO CALCULATE N EXCLUDING BUFFER WITHIN THE NIMBLE MODEL
myHabitat$habitat.rWthBuff <- mask(myHabitat$habitat.r, myHabitat$habitat.poly)
myHabitat$habitat.mxWthBuff <- as.matrix(myHabitat$habitat.rWthBuff)
myHabitat$habitat.mxWthBuff[is.na(myHabitat$habitat.mxWthBuff)] <- 0


### ==== 1.2. GENERATE DETECTORS ====
myDetectors <-  MakeSearchGrid( data = myStudyArea,
                                resolution = myVars$DETECTORS$resolution, 
                                center = TRUE,
                                plot = FALSE)                                         


### ==== 1.3. SCALE DETECTORS & HABITAT & UPPER/LOWER COORDINATES FOR HABITAT WINDOWS ====
scaled <- UTMToGrid( data.sp = myDetectors$detector.sp,
                     grid.sp = myHabitat$habitat.sp,
                     plot.check = FALSE)


ScaledLowerCoords <- UTMToGrid(data.sp = SpatialPoints(lowerHabCoords),
                               grid.sp = myHabitat$habitat.sp)$data.scaled.xy
ScaledUpperCoords <- UTMToGrid(data.sp = SpatialPoints(upperHabCoords),
                               grid.sp = myHabitat$habitat.sp)$data.scaled.xy
ScaledUpperCoords[,2] <- ScaledUpperCoords[,2]+1
ScaledLowerCoords[,2] <- ScaledLowerCoords[,2]-1


### ==== 1.4. CREATE CACHED OBJECTS FOR LESS RESTRICTION ====
DetectorIndexLESS <- GetDetectorIndexLESS( habitat.mx = myHabitat$habitat.mx,
                                           detectors.xy = scaled$data.scaled.xy,
                                           maxDist = 10/res(myHabitat$habitat.r)[1],##--this value should be large enough so that all detections fall withing the restricted window - [EM] what's large enough?
                                           ResizeFactor = 1,
                                           plot.check = FALSE)



## ----------------------------------------------------------------------------------------------
## ------ 2. SIMULATE INPUT DATA  -----
## ----------------------------------------------------------------------------------------------

### ==== 2.1. DETECTOR COVARIATE ====
## simulate autocorrelated rasters
MyDensityRaster <- SimulateAutocorrelatedRaster( sp = myDetectors$detector.sp,
                                                 NRaster = 1,
                                                 phi = myVars$DETECTIONS$phi, ## .001, 1, 1000
                                                 scaled = TRUE,
                                                 focal = NULL,
                                                 plot = FALSE)


##-- fill up raster
## only two categories of detectors (det.prop[1] = 1-det.prop[2])
temp <- MyDensityRaster[]

## Identify the cutoff value between the two categories so that "det.prop" % are
## in category one and 1-det.prop % are in category 2
cutoff <- quantile(temp, myVars$DETECTIONS$det.prop1)

temp[MyDensityRaster[] <= cutoff] <- 1
temp[MyDensityRaster[] > cutoff] <- 2
MyDensityRaster[] <- temp

##-- extract detector covs
myDetectors$detector.sp$detCov <- data.frame(raster::extract(MyDensityRaster, myDetectors$detector.sp))[1]

##--change name to detCov to be consistent over different runs
colnames(myDetectors$detector.sp$detCov) <- "detCov"

##--detector-level covariates
## TURN CATEGORICAL COV INTO MODEL MATRIX
temp <- formula("~detCov - 1")

if(length(unique(unlist(myDetectors$detector.sp$detCov)))>1){
  dummy.cov <- as.data.frame((model.matrix(temp, data.frame(detCov = factor(unlist(myDetectors$detector.sp$detCov))))[,-1]))
  dimnames(dummy.cov)[[2]] <- paste("cov", 1:dim(dummy.cov)[2], sep = "")
  
}else{
  dummy.cov<-matrix(0,nrow = length(myDetectors$detector.sp),ncol = 1)
  dimnames(dummy.cov)[[2]]<-list("cov1")
  dummy.cov<-as.data.frame(dummy.cov)
}



### ==== 2.2. SIMULATE INDIVIDUAL AC LOCATIONS ====
## Sample random (uniform) activity center locations
simulated.ACS <- spsample(x = raster::aggregate(rasterToPolygons(myHabitat$habitat.r,function(x)x>0)),
                          n = myVars$POPULATION$N,
                          type = "random")
simulated.ACS$id <- 1:length(simulated.ACS)


### ==== 2.3. SIMULATE DETECTION ARRAY : y.ar ====

##-- EXPORT EXAMPLE STATE SPACE AND DETECTION FIGURES
detectionSimulationOut <- SimulateDetection(p0 = myVars$DETECTIONS$p0, 
                                            sigma = myVars$DETECTIONS$sigma, 
                                            AC.sp = simulated.ACS, 
                                            detector.sp = myDetectors$detector.sp, 
                                            type = "Bernoulli", 
                                            n.samples = NULL, 
                                            det.covs.p0 = dummy.cov,  ## detector-level covs; active detectors are the intercept
                                            det.betas.p0 = myVars$DETECTIONS$det.betas.p0, 
                                            plot = FALSE)


### ==== 2.4. SUBSET TO INDIVIDUALS DETECTED AT LEAST ONE YEAR/DATA SET ==== 
y <- detectionSimulationOut$y

## CHECK THE NUMBER OF INIDVIDUALS DETECTED
detected <- apply(y,1, function(x) sum(x)>0)

### ==== 2.5. AUGMENT DATA ====

##--augmentation factor that ensures that the total number of individuals 
# (alive + available) is always the same, regardless of the simulation
this.AF <- myVars$POPULATION$M/sum(detected)-1

##--check it: 
sum(detected) * (1+ this.AF) ==  myVars$POPULATION$M

y <- MakeAugmentation(y = y, aug.factor = this.AF, replace.value = 0)
z <- MakeAugmentation( y = rep(1, nrow(detectionSimulationOut$y)),
                       aug.factor = this.AF, replace.value = NA)


### ==== 2.6. TRANSFORM Y TO SPARSE MATRICES =====

SparseY <- GetSparseY(y)


### ==== 2.7.   SET THE INPUT FOR NIMBLE ====
### ==== 2.7.1. Define the nimble model  ====

modelCode <- nimbleCode({
  
  ##---- SPATIAL PROCESS 
  for(i in 1:n.individuals){
    sxy[i,1:2] ~ dbinomPPSingle( lowerHabCoords[1:numHabWindows,1:2],
                                 upperHabCoords[1:numHabWindows,1:2],
                                 mu[1:numHabWindows], 1, numHabWindows)
  }#i
  
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  
  for (i in 1:n.individuals){
    z[i] ~ dbern(psi)
  }#i
  
  ##---- DETECTION PROCESS
  sigma ~ dunif(0, 50)
  alpha <- -1 / (2 * sigma * sigma)
  p0 ~ dunif(0, 1)
  
  
  for (i in 1:n.individuals){
    y[i, 1:nMaxDetectors] ~ dbin_LESSCachedAllSparse(pZero = p0, sxy = sxy[i,1:2], sigma = sigma,
                                                     nbDetections[i], yDets = yDets[i,1:nMaxDetectors],
                                                     detector.xy =  detector.xy[1:n.detectors,1:2],
                                                     trials = trials[1:n.detectors],
                                                     detectorIndex = detectorIndex[1:n.cells,1:maxNBDets], 
                                                     nDetectorsLESS = nDetectorsLESS[1:n.cells],  ResizeFactor = ResizeFactor,
                                                     maxNBDets = maxNBDets, habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],maxDist=MaxDist,
                                                     indicator = z[i])
  }#i
  
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
  
})


### ==== 2.7.2. Define the constants to be used in the model ====
# Set the list of model constants
nimConstants <- list(
  x.max = habQualityDims[1],
  y.max =  habQualityDims[2],
  n.detectors = dim(y)[2], 
  n.individuals = dim(y)[1],
  numHabWindows = nHabCells,
  nMaxDetectors = SparseY$nMaxDetectors,
  y.maxDet = dim(DetectorIndexLESS$habitatID)[1],
  x.maxDet = dim(DetectorIndexLESS$habitatID)[2],
  ResizeFactor = DetectorIndexLESS$ResizeFactor,
  n.cellsSparse = dim(DetectorIndexLESS$detectorIndex)[1],
  maxNBDets = DetectorIndexLESS$maxNBDets,
  n.cells = dim(DetectorIndexLESS$detectorIndex)[1],
  MaxDist = myVars$DETECTIONS$sigma *  5 
)


### ==== 2.7.3. Define the data to be used in the model ====

# Set the list of data
nimData <- list(
  lowerHabCoords = ScaledLowerCoords,
  upperHabCoords = ScaledUpperCoords,
  z = z,
  y = SparseY$y[,,1],
  detector.xy = scaled$data.scaled.xy,
  trials = rep(1, nrow(scaled$data.scaled.xy)),  
  mu = rep(1, nHabCells),
  yDets = SparseY$yDets[,,1],
  nbDetections = SparseY$nbDetections[,1],
  detectorIndex = DetectorIndexLESS$detectorIndex,
  nDetectorsLESS = DetectorIndexLESS$nDetectorsLESS,
  habitatIDDet = DetectorIndexLESS$habitatID
)


### ==== 2.7.4. Define the initial values to be used in the model  ====
## Initialise AC locations
sxy <- MakeInitXY( y = y,
                   habitat.mx = myHabitat$habitat.mx,
                   detector.xy = scaled$data.scaled.xy,
                   IDCells.mx = myHabitat$IDCells.mx,
                   grid.xy = scaled$grid.scaled.xy)

## Initialise z values
z.init <- ifelse(!is.na(z), NA, 1)
z.init[!is.na(z.init)] <- rbinom(sum(!is.na(z.init)), size = 1, prob = 0.5)


nimInits <- list(sigma = runif(1,0,10),
                 p0 = runif(1,0,1),   
                 z = z.init,           
                 psi = runif(1,0,1),
                 sxy = sxy[ , ,1]      
)


### ==== 2.7.5. Define the parameters to be monitored ====
nimParams <- c("N","sxy", "z", "sigma","p0", "psi")



## ----------------------------------------------------------------------------------------------
## ------ 3. SCR - A TEST RUN
## ----------------------------------------------------------------------------------------------

### ==== 3.2. CREATE & COMPILE THE NIMBLE MODEL OBJECT & MCMC CONFIGURATION ====

model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = FALSE,       
                      calculate = FALSE)  

cmodel <- compileNimble(model)


MCMCconf <- configureMCMC(model = model,
                          monitors = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1) 

MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)


### ==== 3.3.RUN THE MCMC ====
myNimbleOutput <- runMCMC(mcmc = cMCMC,
                          nburnin = 1, 
                          niter = 100,  
                          nchains = 3, 
                          samplesAsCodaMCMC = TRUE)

