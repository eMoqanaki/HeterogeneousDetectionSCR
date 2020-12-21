#' @title Function to simulate individual detections within an SCR framework, under the influence of 
#' different individual and spatial (detector-level) covariates.
#'
#' @description
#' \code{SimulateDetection} returns a list object with \code{} 
#' 
#' @param p0 \code{Numeric} variable denoting the intercept value of the half-normal detection function describing the 
#' decay of detection probability with increasing distance from the AC.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the half-normal detection function.
#' @param AC.sp A \code{SpatialPointsDataFrame} object with individual activity center locations.
#' @param detector.sp  A \code{SpatialPointsDataFrame} object with detectors' locations.
#' @param id.covs.p0 A \code{DataFrame} with the individual-level covariates to affect p0.
#' @param id.betas.p0 A \code{vector} object with the coefficients associated to each individual covariate. 
#' @param det.covs.p0 A \code{DataFrame} with the detector-level covariates  to affect p0.
#' @param det.betas.p0 A \code{vector} object with the coefficients associated to each detector covariate. 
#' @param id.covs.sigma A \code{DataFrame} with the individual-level covariates to affect sigma. 
#' @param id.betas.sigma A \code{vector} object with the coefficients associated to each individual covariate. 
#' @param link.sigma A \code{character} to specify the link-function to be used when modelling the effect of covariates on sigma (default is exponential to avoid negative sigmas)
#' @param n.samples A \code{numeric} object denoting the total number of samples to be returned.
#' @param alive.ids A \code{vector} vector denoting the ids of the individuals that are alive. They should match the AC.sp. 
#' @param alpha A \code{probability} denoting the identification probability (implementation of partial identification for Poisson simulations)
#' @param type A \code{character} to specify the type of detections to be modelled; can be binary (\code{"Bernoulli"} or \code{"Binomial})
#' or it can be derived from a count process (\code{"Poisson"} or \code{"BinomFomPoisson}).
#' @param plot.check A \code{logical} for whether (\code{TRUE}) or not (\code{FALSE}) plots are to be generated during simulations.
#' @param seed A seed number for results to be reproducible (default is NULL)

#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' det.sim<-SimulateDetection(p0=0.6,sigma=1200,AC.sp=AC.sp,detector.sp=detector.sp,covs=covs, betas=betas,plot=TRUE)

SimulateDetection <- function( p0 = NULL 
                             , a0 = NULL
                             , sigma
                             , AC.sp
                             , detector.sp
                             , type = "Bernoulli"
                             
                                ## EFFECTS ON p0
                             # Detector Fixed Effects
                             , det.covs.p0 = NULL
                             , det.betas.p0 = NULL
                             # Detector Random Effects
                             , det.rand.sd.p0 = NULL
                             # Detector Offset Effects
                             , det.offset.p0 = NULL
                             
                             # Individual Fixed Effects
                             , id.covs.p0 = NULL
                             , id.betas.p0 = NULL
                             # Individual Random Effects
                             , id.rand.sd.p0 = NULL
                                  
                                ## EFFECTS ON SIGMA
                             , link.sigma = "exponential"
                             # Individual Fixed Effects
                             , id.covs.sigma = NULL
                             , id.betas.sigma = NULL
                             # Individual Random Effects
                             , id.rand.sd.sigma= NULL
                             
                                ## OTHERS
                             , n.samples = NULL
                             , alive.ids = NULL
                             , alpha = NULL
                             , plot.check = TRUE
                             , habitat.poly = NULL
                             , seed = NULL
                             , typeD = 1)
   {
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== I.CLEAN & SET-UP THE INPUT DATA ====
   projection(AC.sp) <- projection(detector.sp)
   
   if(!is.null(p0)){myP0 <- p0}                                  ## For plotting purpose
   if(!is.null(a0)){myP0 <- 1-exp(-(a0/(2*pi*sigma^2)))}         ## For plotting purpose
   mySIGMA <- sigma                                              ## For plotting purpose
   
   #--- Create separate detectors and sub-detectors SpatialPointsDataFrame when needed
   if(type == "Binomial" | type == "Bernoulli"){
      #--- Subset detector.sp to main detectors only
      sub.detector.sp <- detector.sp
      temp <- unique(detector.sp@data[ ,c("main.cell.id","main.cell.x","main.cell.y")])
      temp <- temp[order(temp$main.cell.id), ]
      detector.sp <- SpatialPointsDataFrame(temp[ ,c("main.cell.x","main.cell.y")], data = data.frame(temp), proj4string = CRS(projection(detector.sp)))
      
      #--- Extract detector-specific number of trials
      n.trials <- unlist(lapply(detector.sp@data$main.cell.id,function(x){sum(sub.detector.sp@data$main.cell.id==x)}))
      n.trials <- matrix(n.trials, length(AC.sp), length(detector.sp), byrow = TRUE) 
      if(type == "Bernoulli"){n.trials[] <- 1}
      }#if
   
   #--- Check if detectors and covariates dimensions match
   if(!is.null(det.covs.p0)){
      if(dim(det.covs.p0)[1] != length(detector.sp)){
         stop("Detectors' covariates dataframe dimensions do not match detectors SpatialPointDataFrame dimensions!")                
         }#if
      }#if
   if(!is.null(det.offset.p0)){
      if(length(det.offset) != length(detector.sp)){
         stop("Detectors' offset covariate vector length does not match detectors SpatialPointDataFrame dimensions!")                
         }#if
      }#if
   if(!is.null(id.covs.p0)){
      if(dim(id.covs.p0)[1] != length(AC.sp)){
         stop("Individual-level covariates for p0 do not match Individual ACs SpatialPointDataFrame dimensions!")                
      }#if
   }#if
   if(!is.null(id.covs.sigma)){
      if(dim(id.covs.sigma)[1] != length(AC.sp)){
         stop("Individual-level covariates for sigma do not match Individual ACs SpatialPointDataFrame dimensions!")                
      }#if
   }#if
   
   #--- Check if arguments are not contradictory
   if(!is.null(p0) & !is.null(a0)){
      stop("Unauthorized parameters : you cannot provide both p0 and a0, pick one!")  
      }#if
   
   #--- Occasionally rownames pose a problem; fix now:
   dimnames(detector.sp@data)[[1]] <- dimnames(detector.sp@coords)[[1]] <- 1:length(detector.sp)
   dimnames(AC.sp@data)[[1]] <- dimnames(AC.sp@coords)[[1]] <- 1:length(AC.sp)
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== II.CALCULATE INDIVIDUAL-SPECIFIC sigmas ====
      ### ---- 1.Individual-specific covariates effects on sigma ----
   # Fixed effects
   id.fixed.effects.sigma <- 0
   if(!is.null(id.covs.sigma)){
      temp <- formula(paste("~", paste(names(id.covs.sigma),collapse="+"),"-1",sep=""))
      Xmat <- model.matrix(temp, id.covs.sigma)
      id.fixed.effects.sigma <- Xmat%*%id.betas.sigma
      id.fixed.effects.sigma <- do.call(cbind, lapply(1:length(detector.sp), function(x)id.fixed.effects.sigma))
      X.ind.p0 <- data.frame(Xmat)
   }#if
   
   # Random effects
   id.random.effects.sigma <- 0
   if(!is.null(id.rand.sd.sigma)){
      id.random.effects.sigma <- rnorm(n = length(AC.sp), mean = 0, sd = id.rand.sd.sigma)
      id.random.effects.sigma <- do.call(cbind, lapply(1:length(detector.sp), function(x)id.random.effects.sigma))
      }#if
   
      ### ---- 2.Calculate individual-specific sigma ----
   if(link.sigma=="exponential"){
      intercept.sigma <- rep(log(sigma),length(AC.sp))
      sigma <- exp(intercept.sigma + id.fixed.effects.sigma + id.random.effects.sigma)
   }else{
      intercept.sigma <- rep(sigma,length(AC.sp))
      sigma <- intercept.sigma + id.fixed.effects.sigma + id.random.effects.sigma
   }
   sigma <- matrix(sigma, length(AC.sp), length(detector.sp), byrow = FALSE)
   individual.sigma <- sigma[,1]
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== III.CALCULATE DETECTOR & INDIVIDUAL-SPECIFIC p0s ====
      ### ---- 1.Detector-specific covariates effects on p0 ----
   # Fixed effects
   det.fixed.effects.p0 <- 0
   if(!is.null(det.covs.p0)){
      temp <- formula(paste("~", paste(names(det.covs.p0), collapse = "+"), "-1", sep = ""))
      Xmat <- model.matrix(temp, det.covs.p0)
      det.fixed.effects.p0 <- Xmat%*%det.betas.p0
      det.fixed.effects.p0 <- do.call(rbind, lapply(1:length(AC.sp), function(x)t(det.fixed.effects.p0)))
      X.det.p0 <- data.frame(Xmat)
      }#if
   
   # Random effects
   det.random.effects.p0 <- 0
   if(!is.null(det.rand.sd.p0)){
      det.random.effects.p0 <- rnorm(n = length(detector.sp), mean = 0, sd = det.rand.sd.p0)
      det.random.effects.p0 <- do.call(rbind, lapply(1:length(AC.sp), function(x)t(det.random.effects.p0)))
      }#if
   
      ### ---- 2.Individual-specific covariates effects on p0 ----
   # Fixed effects
   id.fixed.effects.p0 <- 0
   if(!is.null(id.covs.p0)){
      temp <- formula(paste("~", paste(names(id.covs.p0),collapse="+"),"-1",sep=""))
      Xmat <- model.matrix(temp, id.covs.p0)
      id.fixed.effects.p0 <- Xmat%*%id.betas.p0
      id.fixed.effects.p0 <- do.call(cbind, lapply(1:length(detector.sp), function(x)id.fixed.effects.p0))
      X.ind.p0 <- data.frame(Xmat)
      }#if
   
   # Random effects
   id.random.effects.p0 <- 0
   if(!is.null(id.rand.sd.p0)){
      id.random.effects.p0 <- rnorm(n = length(AC.sp), mean = 0, sd = id.rand.sd.p0)
      id.random.effects.p0 <- do.call(cbind, lapply(1:length(detector.sp), function(x)id.random.effects.p0))
      }#if
   
   inv.logit <- function(inValues) {
      1.0 / (1.0 + exp(-inValues))
   }
   
      ### ---- 3.Calculate individual & detector-specific p0 ----
   if(!is.null(a0)){p0 <- 1-exp(-(a0/(2*pi*sigma^2)))}else{matrix(p0,length(AC.sp),length(detector.sp))} 
   
   if(type == "Bernoulli" | type == "Binomial"){
      intercept.p0 <- matrix(logit(p0), length(AC.sp), length(detector.sp))
      p0 <- inv.logit(intercept.p0 + id.fixed.effects.p0 + id.random.effects.p0 + det.fixed.effects.p0 + det.random.effects.p0)
      }#if
   if(type == "Poisson"){
      intercept.p0 <- matrix(log(p0), length(AC.sp), length(detector.sp))
      p0 <- exp(intercept.p0 + id.fixed.effects.p0 + id.random.effects.p0 + det.fixed.effects.p0 + det.random.effects.p0)
      }#if

      ### ---- 4.Additional offset effect on p0 ----
   if(!is.null(det.offset.p0)){
      det.offset <- do.call(rbind, lapply(1:length(AC.sp), function(x)t(det.offset)))
      p0 <- 1-(1-p0)^det.offset
      }#if
   
   individual.p0 <- p0
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== IV.CALCULATE INDIVIDUAL & DETECTOR-SPECIFIC p ====
   #--- Calculate distance matrix between detectors and ACs
   if(typeD == 1){D <- gDistance(detector.sp, AC.sp, byid=TRUE)}
   if(typeD == 2){  
      x1 <- detector.sp$main.cell.x
      y1 <- detector.sp$main.cell.y
      x2 <- AC.sp$x
      y2 <- AC.sp$y
      D <- matrix(unlist(lapply(1:length(x2), function(dd){sqrt((x1 - x2[dd])^2 + (y1 - y2[dd])^2) })), length(x2), length(x1), byrow = TRUE)
      }#if  
 
   #--- Calculate Individual & Detector-specific detection probability
   P <- p0*exp(-D*D/(2*sigma*sigma))
   
   #--- Set dead individuals detection probability to 0
   if(!is.null(alive.ids)){
      P2 <- matrix(0, dim(P)[1], dim(P)[2])
      P2[alive.ids, ] <- P[alive.ids, ]
      P <- P2
      }
   
   ##-------------------------------------------------------------------------------------------------------------------
   ## ==== V.BERNOULLI or BINOMIAL DETECTION PROCESS ====
   if(type == "Binomial" | type == "Bernoulli"){      
      #-- Set the seed for replication
      y.seed <- P
      if(!is.null(seed)){set.seed(seed)}
      y.seed[] <- sample(1000000, length(n.trials), replace = TRUE)
      #-- Sample number of individual detections per main detector
      temp <- abind(n.trials, P, y.seed, along = 3)
      y <- apply(temp, c(1,2), function(x){
         set.seed(x[3])
         rbinom(1, x[1], x[2])})
         }#if              
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== VI.POISSON DETECTION PROCESS ====
   if(type == "Poisson"){
      #-- Set the seed for replication
      y.seed <- P
      if(!is.null(seed)){set.seed(seed)}
      y.seed[] <- sample(1000000, length(P), replace = TRUE)
      
      #-- Sample number of samples found per detector
      temp <- abind(y.seed, P,along = 3)
      y <- apply(temp, c(1,2), function(x){
         set.seed(x[1])
         rpois(1, x[2])}) 
      
      #--- Resample to match n.samples per individual
      if(!is.null(n.samples)){
         if(n.samples < sum(y)){
            y.sub <- y       
            y.sub[] <- apply(rmultinom(n.samples, 1, p = y/sum(y)), 1, sum)
            p0 <- p0*sum(y.sub)/sum(y)                                        # to output correct parameter after subsampling
            y <- y.sub
            }#if
         }#if
      
      #--- Partial identification: sample individuals unidentified
      ind.y <- y
      if(!is.null(alpha)){
         ind.y[] <- unlist(lapply(y, function(x){rbinom(1, x, alpha)}))
         count.y <- apply(y, 2, sum)                                         # Counts data: known and unknown ids
         y <- ind.y                                                          # Identified IDs only
         }#if
      }#if
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== VII.BINOMIAL FROM POISSON DETECTION PROCESS ====
   if(type == "BinomFromPoisson"){ 
      #-- Set the seed for replication
      y.seed <- P
      if(!is.null(seed)){set.seed(seed)}
      y.seed[] <- sample(1000000, length(P), replace = TRUE)
         
      #-- Sample number of samples found per detector
      temp <- abind(y.seed, P,along = 3)
      y <- apply(temp, c(1,2), function(x){
               set.seed(x[1])
               rpois(1, x[2])}) 
         
      #--- Resample to match n.samples per individual
      if(!is.null(n.samples)){
         if(n.samples < sum(y)){
            y.sub <- y       
            y.sub[] <- apply(rmultinom(n.samples, 1, p = y/sum(y)), 1, sum)
            p0 <- p0*sum(y.sub)/sum(y)                                        # to output correct parameter after subsampling
            y <- y.sub
            }
         }
         
      #--- Partial identification: sample individuals unidentified
      ind.y <- y
      if(!is.null(alpha)){                              
         ind.y[] <- unlist(lapply(y, function(x){rbinom(1, x, alpha)}))
         count.y <- apply(y, 2, sum)                                         # Counts data: known and unknown ids
         y <- ind.y                                                          # Identified IDs only
         }
   
      #--- Aggregate Poisson counts to Bernoulli detections at the sub-detector level
      y.binom <- y
      y.binom[y>0] <- 1
      
      #--- Aggregate to main detectors
      y.binom <- t(aggregate(t(y.binom), by=list(detector.sp@data$main.cell.id), FUN=sum))
      row.names(y.binom) <- 1:dim(y.binom)[1]
      y.binom <- y.binom[-c(1), ]
      rownames(y.binom) <- 1:dim(y.binom)[1]
      y <- y.binom
      
      #--- Subset detector.sp to main detectors only
      sub.detector.sp <- detector.sp
      temp <- unique(detector.sp@data[ ,c("main.cell.id","main.cell.x","main.cell.y")])
      temp <- temp[order(temp$main.cell.id), ]
      detector.sp <- SpatialPointsDataFrame(temp[ ,c("main.cell.x","main.cell.y")], data=data.frame(temp), proj4string=CRS(projection(detector.sp)))
      }#if
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== VIII.SPATIAL DEAD RECOVERY PROCESS ====
   if(type == "Recovery"){ 
      addedP <- rep(1, dim(P)[1])
      if(length(alive.ids) > 0){addedP[alive.ids] <- 0}
      P <- cbind(addedP,P)
      y <- apply(P, 1, function(x){which(rmultinom(1,1,x)==1)})
      y <- (y-1)
      }#if 
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== IX.OUTPUT LIST ====
   #--- Fix the names
   if(!is.vector(y)){dimnames(y) <- list(1:dim(y)[1], 1:dim(y)[2])}  #---DO NOT REMOVE! :)
   dimnames(AC.sp@coords) <- list(c(1:length(AC.sp)), c("x","y"))
    
   #--- Slim to individuals detected
   y.all <- y
   if(length(dim(y)) > 1){
      detected <- apply(y, 1, max) > 0
      if(!any(detected)){detected <- rep(TRUE, dim(y)[1])}
      y <- y[detected, ]
      }
   if(length(dim(y)) == 0){
      detected <- which(y < dim(P)[2])
      if(!any(detected)){detected <- rep(TRUE, length(y))}
      y <- y[detected]
      }
   
   #--- Output list
   out <- list( y = y        
              , y.all = y.all
              , D = D
              , p0 = p0
              , sigma = sigma
              , id.covs.p0 = data.frame(id.covs.p0[detected, ])
              , id.betas.p0 = id.betas.p0
              , det.covs.p0 = det.covs.p0
              , det.betas.p0 = det.betas.p0
              , id.covs.sigma = data.frame(id.covs.sigma[detected, ])
              , id.betas.sigma = id.betas.sigma
              , n.samples = n.samples
              , alpha = alpha
              , individual.sigma = individual.sigma
              , individual.p0 = individual.p0)  
   
   if(type == "Poisson" & !is.null(alpha)){out$count.y <- count.y}
   if(type == "Binomial" |type == "Bernoulli" ){out$n.trials <- n.trials[1, ]}
   
   ##------------------------------------------------------------------------------------------------------------------- 
   ## ==== X.PLOT DETECTIONS ====
   if(plot.check){
      ### ---- 1.Plot individual cov effects on p0 ----
      if(!is.null(id.covs.p0) & type != "BinomFromPoisson" & type != "Recovery"){
         
         AC.sp@data <- id.covs.p0
         
         for(i in 1:dim(id.covs.p0)[2]){
            myCov <- id.covs.p0[ ,i]
            myValues <- unique(myCov) 
            myValues <- myValues[order(myValues)]
            myCols <- rev(heat.colors(n = length(myValues), alpha = 1))
            
            par(mfrow=c(1,2))
            
            ## Plot the effect of covariate i on p0
            myBeta <- id.betas.p0[i]
            myX <- seq(min(myValues), max(myValues), length.out = 1000)
            
            if(type == "Poisson"){
               myY <- exp(log(myP0) + myBeta*myX)
               plot(myX, myY, xlab = names(id.covs.p0[i]), ylab = "lambda0", ylim = c(0,max(myY)), pch=19, col = rev(heat.colors(length(myX))))
               if(!is.null(id.rand.sd.p0)){
                  myY.up <- exp(log(myP0 + qnorm(0.025,0,id.rand.sd.p0)) + myBeta*myX)
                  myY.down <- exp(log(myP0 + qnorm(0.975,0,id.rand.sd.p0)) + myBeta*myX)
                  points(myX, myY.up,  pch=19, cex = 0.2, col = rev(heat.colors(length(myX))))
                  points(myX, myY.down,  pch=19, cex = 0.2, col = rev(heat.colors(length(myX))))
                  }#if
               }#if
            
            if(type != "Poisson"){
               myY <- inv.logit(logit(myP0) + myBeta*myX)
               plot(myX, myY, xlab = names(id.covs.p0[i]), ylab = "p0", ylim = c(0,1),  pch=19, col = rev(heat.colors(length(myX))))
               if(!is.null(id.rand.sd.p0)){
                  myY.up <- inv.logit(logit(myP0 + qnorm(0.025,0,id.rand.sd.p0)) + myBeta*myX)
                  myY.down <- inv.logit(logit(myP0 + qnorm(0.975,0,id.rand.sd.p0)) + myBeta*myX)
                  points(myX, myY.up,  pch=19, cex = 0.5, col = rev(heat.colors(length(myX))))
                  points(myX, myY.down,  pch=19, cex = 0.5, col = rev(heat.colors(length(myX))))
                  }#if
               }#if
            
            ## Plot the Detections
            plot(AC.sp, main = names(id.covs.p0)[i])
            if(!is.null(habitat.poly)){plot(habitat.poly, col = rgb(t(col2rgb("forestgreen")/255),alpha = 0.6), add = TRUE)}
            plot(detector.sp, col = "gray40", pch = 19, cex=0.6, add=TRUE)
            plot(AC.sp, add = TRUE)
            
            for(j in 1:length(myValues)){
               myIds <- which(AC.sp@data[ ,i] == myValues[j])
               
               lapply(myIds,function(x){
                  this.row <- y.all[x, ]
                  if(sum(this.row)>0){
                     plot(AC.sp[x, ], col = myCols[j], pch = 19, add = TRUE)
                     this.det <- detector.sp[this.row>0, ]
                     segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                             , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                             , col = myCols[j])
                     }
                  })
               }
            }
         }
      
      ### ---- 2.Plot detector cov effects on p0 ----
      if(!is.null(det.covs.p0) & type != "BinomFromPoisson" & type != "Recovery"){
         
         detector.sp@data <- det.covs.p0
         
         for(i in 1:dim(det.covs.p0)[2]){
            myCov <- det.covs.p0[ ,i]
            myValues <- unique(myCov) 
            myValues <- myValues[order(myValues)]
            myCols <- rev(terrain.colors(n = length(myValues), alpha = 1))
            
            par(mfrow=c(1,2))
            
            ## Plot the effect of covariate i on p0
            myBeta <- det.betas.p0[i]
            myX <- seq(min(myValues),max(myValues), length.out = 1000)
            if(type == "Poisson"){
               myY <- exp(log(myP0) + myBeta*myX)
               plot(myX, myY, xlab = names(det.covs.p0[i]), ylab = "lambda0", ylim = c(0,max(myY)), pch=19, col = rev(terrain.colors(length(myX))))
               if(!is.null(det.rand.sd.p0)){
                  myY.up <- exp(log(myP0 + qnorm(0.025,0,det.rand.sd.p0)) + myBeta*myX)
                  myY.down <- exp(log(myP0 + qnorm(0.975,0,det.rand.sd.p0)) + myBeta*myX)
                  points(myX, myY.up,  pch=19, cex = 0.2, col =rev(terrain.colors(length(myX))))
                  points(myX, myY.down,  pch=19, cex = 0.2, col = rev(terrain.colors(length(myX))))
                  }#if
               }#if
            if(type != "Poisson"){
               myY <- inv.logit(logit(myP0) + myBeta*myX)
               plot(myX, myY, xlab = names(det.covs.p0)[i], ylab = "p0", ylim = c(0,1), pch=19, col = rev(terrain.colors(length(myX))))
               if(!is.null(det.rand.sd.p0)){
                  myY.up <- inv.logit(logit(myP0 + qnorm(0.025,0,det.rand.sd.p0)) + myBeta*myX)
                  myY.down <- inv.logit(logit(myP0 + qnorm(0.975,0,det.rand.sd.p0)) + myBeta*myX)
                  points(myX, myY.up,  pch=19, cex = 0.5, col = rev(terrain.colors(length(myX))))
                  points(myX, myY.down,  pch=19, cex = 0.5, col = rev(terrain.colors(length(myX))))
                  }#if
               }#if
            
            ## Plot the Detections
            plot(AC.sp, main = names(det.covs.p0)[i])
            if(!is.null(habitat.poly)){plot(habitat.poly, col = rgb(t(col2rgb("forestgreen")/255),alpha = 0.6), add = TRUE)}
            
            for(j in 1:length(myValues)){
               myIds <- which(detector.sp@data[ ,i] == myValues[j])
               plot(detector.sp[myIds,], col = myCols[j], pch=19, cex=2, add=TRUE)
               }#j
            
            plot(AC.sp, add = TRUE)
            plot(AC.sp[detected, ], col = "red", pch = 19, add = TRUE)
            
            lapply(which(detected), function(x){
               this.row <- y.all[x, ]
               this.det <- detector.sp[this.row>0, ]
               if(length(this.det)>0){
                  segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                            , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                            , col = "red")
               }})
            }#i
         }#if
      
      ### ---- 3.Plot individual cov effects on sigma -----
      if(!is.null(id.covs.sigma) & type != "BinomFromPoisson" & type != "Recovery"){
         
         AC.sp@data <- id.covs.sigma
         
         for(i in 1:dim(id.covs.sigma)[2]){
            myCov <- id.covs.sigma[ ,i]
            myValues <- unique(myCov) 
            myValues <- myValues[order(myValues)]
            myCols <- rev(heat.colors(n = length(myValues), alpha = 1))
            
            par(mfrow = c(1,2))
            ## Plot the effect of covariate i on sigma
            myBeta <- id.betas.sigma[i]
            myX <- seq(min(myValues), max(myValues), length.out = 1000)
            
            if(link.sigma == "exponential"){
               myY <- exp(log(mySIGMA) + myBeta*myX)
               plot(myX, myY, xlab = names(id.covs.sigma[i]), ylab = "sigma", ylim = c(0,max(myY)), pch=19, col = rev(heat.colors(length(myX))))
               if(!is.null(id.rand.sd.sigma)){
                  myY.up <- exp(log(mySIGMA + qnorm(0.025,0,id.rand.sd.sigma)) + myBeta*myX)
                  myY.down <- exp(log(mySIGMA + qnorm(0.975,0,id.rand.sd.sigma)) + myBeta*myX)
                  points(myX, myY.up, cex = 0.5, pch=19, col = rev(heat.colors(length(myX))))
                  points(myX, myY.down, cex = 0.5, pch=19, col = rev(heat.colors(length(myX))))
                  }#if
               }#if
            
            if(link.sigma != "exponential"){
               myY <- mySIGMA + myBeta*myX
               plot(myX, myY, xlab = names(id.covs.sigma[i]), ylab = "sigma", ylim = c(0,max(myY)), pch=19, col = rev(heat.colors(length(myX))))
               if(!is.null(id.rand.sd.sigma)){
                  myY.up <- mySIGMA + qnorm(0.025,0,id.rand.sd.sigma) + myBeta*myX
                  myY.down <- mySIGMA + qnorm(0.975,0,id.rand.sd.sigma) + myBeta*myX
                  points(myX, myY.up, cex = 0.5, pch=19, col = rev(heat.colors(length(myX))))
                  points(myX, myY.down, cex = 0.5, pch=19, col = rev(heat.colors(length(myX))))
                  }#if
               }#if
            


            ## Plot the Detections
            plot(AC.sp, main = names(id.covs.sigma)[i])
            if(!is.null(habitat.poly)){plot(habitat.poly, col = rgb(t(col2rgb("forestgreen")/255),alpha = 0.6), add = TRUE)}
            plot(detector.sp, col="gray40", pch=19, cex=0.6, add=TRUE)
            plot(AC.sp, add = TRUE)
            
            for(j in 1:length(myValues)){
               myIds <- which(AC.sp@data[ ,i] == myValues[j])
               
               lapply(myIds,function(x){
                  this.row <- y.all[x, ]
                  this.det <- detector.sp[this.row>0, ]
                  if(length(this.det)>0){
                     plot(AC.sp[x, ], col = myCols[j], pch = 19, add = TRUE)
                     segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                               , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                               , col = myCols[j])}})
            }
         }
      }
      
      ### ---- 4.Plot compensatory effect of sigma on p0 ----
      if(!is.null(a0) & type != "BinomFromPoisson" & type != "Recovery"){

         par(mfrow=c(1,2))
            
         ## Plot the compensatory effect between sigma & p0
         if(length(unique(individual.sigma)) <= 3){mySigma <- unique(individual.sigma)}
         if(length(unique(individual.sigma)) > 3){mySigma <- as.vector(c(quantile(individual.sigma, 0.025), mean(individual.sigma), quantile(individual.sigma, 0.975)))}
         myp0 <- 1-exp(-(a0/(2*pi*mySigma^2)))
         myCols <- rev(heat.colors(n = length(mySigma), alpha = 1))
         myCols.fill <- rev(heat.colors(n = length(mySigma), alpha = 0.5))
         
         temp.X <- seq(0, 3*max(mySigma), length.out = 1000)
         temp.Y1 <- myp0[1]*exp(-(temp.X^2)/(2*mySigma[1]^2))
         temp.Y2 <- myp0[2]*exp(-(temp.X^2)/(2*mySigma[2]^2))
         temp.Y3 <- myp0[3]*exp(-(temp.X^2)/(2*mySigma[3]^2))
             
         myX <- c(-rev(temp.X), temp.X)
         myY1 <- c(rev(temp.Y1), temp.Y1)                                                                 
         myY2 <- c(rev(temp.Y2), temp.Y2)                                                                 
         myY3 <- c(rev(temp.Y3), temp.Y3) 
         
         plot(myX, myY1, type="n", xlab = "distance", ylab = "p[i,j]", xlim = c(-3*max(mySigma), 3*max(mySigma)), ylim = c(0,1), col = myCols[1])
         polygon(x = c(myX, rev(myX)), y = c(myY1, rep(0,length(myY1))), border = myCols[1], col = myCols.fill[1])        
         polygon(x = c(myX, rev(myX)), y = c(myY2, rep(0,length(myY2))), border = myCols[2], col = myCols.fill[2])        
         polygon(x = c(myX, rev(myX)), y = c(myY3, rep(0,length(myY3))), border = myCols[3], col = myCols.fill[3])        
         
         ## Plot the Detections
         myValues <- unique(individual.sigma) 
         myValues <- myValues[order(myValues)]
         myCols <- rev(heat.colors(n = length(myValues), alpha = 1))
         
         plot(AC.sp, main = "p0/sigma compensation")
         if(!is.null(habitat.poly)){plot(habitat.poly, col = rgb(t(col2rgb("forestgreen")/255),alpha = 0.6), add = TRUE)}
         plot(detector.sp, col = "gray40", pch = 19, cex=0.6, add=TRUE)
         plot(AC.sp, add = TRUE)
         
         for(j in 1:length(myValues)){
            myIds <- which(individual.sigma == myValues[j])
            
            lapply(myIds,function(x){
               this.row <- y.all[x, ]
               if(sum(this.row)>0){
                  plot(AC.sp[x, ], col = myCols[j], pch = 19, add = TRUE)
                  this.det <- detector.sp[this.row>0, ]
                  segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                          , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                          , col = myCols[j])
                  }#if
               })#lapply
            }#j
         }#if
      
      ### ---- 5.Plot individual detections (no covariates) ----
      if(is.null(id.covs.sigma) & is.null(id.covs.p0) & is.null(det.covs.p0) & type != "Recovery" | type == "BinomFromPoisson")
         {
         if(!is.null(det.rand.sd.p0) & !is.null(id.rand.sd.sigma) | !is.null(id.rand.sd.p0) & !is.null(id.rand.sd.sigma)){
            par(mfrow = c(1,3))
            hist(individual.p0, col = "gray40",  main = "Random variation in p0")
            abline(v=myP0, col="red", lwd=2)
            hist(individual.sigma, col = "gray40", main = "Random variation in sigma")
            abline(v=mySIGMA, col="red", lwd=2)
            }#if
         if(!is.null(det.rand.sd.p0) & is.null(id.rand.sd.sigma) | !is.null(id.rand.sd.p0) & is.null(id.rand.sd.sigma)){
            par(mfrow = c(1,2))
            hist(individual.p0, col = "gray40", main = "Random variation in p0")
            abline(v=myP0, col="red", lwd=2)
            }#if
         if(is.null(det.rand.sd.p0) & is.null(id.rand.sd.p0) & !is.null(id.rand.sd.sigma)){
            par(mfrow = c(1,2))
            hist(individual.sigma, col = "gray40", main = "Random variation in sigma")
            abline(v=mySIGMA, col="red", lwd=2)
            }#if
         
         plot(AC.sp)
         if(!is.null(habitat.poly)){plot(habitat.poly, col = rgb(t(col2rgb("forestgreen")/255),alpha = 0.6), add = TRUE)}
         plot(detector.sp, col = "gray40", pch = 19, cex = 0.6, add = TRUE)
         plot(AC.sp, add = TRUE)
         
         plot(AC.sp[detected, ], col = "red", pch = 19, add = TRUE)
         
         lapply(which(detected), function(x){
            this.row <- y.all[x, ]
            this.det <- detector.sp[this.row>0, ]
            if(length(this.det)>0){
               segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                       , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                       , col = "red")
            }})
         }
      
      ### ---- 6.Plot individual dead recoveries ----
      if(type == "Recovery"){
         plot(AC.sp)
         if(!is.null(habitat.poly)){plot(habitat.poly, col = rgb(t(col2rgb("forestgreen")/255),alpha = 0.6), add = TRUE)}
         plot(detector.sp, col="gray40", pch=19, cex=0.6, add = TRUE)
         plot(AC.sp, add = TRUE)
         
         recovered <- which(y.all != 0)
         plot(AC.sp[recovered,], col = "red", pch = 19, add = TRUE)
         
         lapply(recovered, function(x){
            this.det <- detector.sp[y.all[x], ]
            segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                    , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                    , col = "red")})
         }
      }#if
   ##------------------------------------------------------------------------------------------------------------------- 
   return(out)
   }
