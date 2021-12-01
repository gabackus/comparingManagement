rm(list=ls())
require(pracma)
require(extraDistr)
require(parallel)

commSetup <- function(S=64, L=512, W=8,
                      zo=NULL, gam=NULL, sig=NULL, A=NULL, m=1,
                      gamMean=2.5, gamSD=2.5,
                      sigMean=5, sigSD=5,
                      lam=-2.7, B=10,
                      compType="lottery",
                      XW=seq(129,384),
                      temp2d=NULL,
                      tempLow=9.78, tempHigh=30.22,
                      tempRev=F,
                      tempGH=1, tempGSD=0, tempLSD=0,
                      years=1000,
                      tempY=NULL,
                      tau=0.04,
                      tempYAC=0.767, tempYSD=0.1639,
                      Tau=NULL,
                      tempYXGH=1, tempYXGSD=0, tempYXLSD=0,
                      tempYXGLH=1, tempYXGLSD=0,
                      Q=NULL,
                      QMean=8,
                      QGH=1, QGSD=0, QLSD=0){
  
  # This sets up the basic structure of the model.
  # It creates a list of biological parameters in P (S randomized species and their parameters)
  # and it creates a list of environmental parameters in X
  
  # S:         Total number of species created with hetSetup.
  # L:         Total number of patches in the metacommunity.
  # W:         Number of microhabitats in each patch.
  
  # zo:        A vector of pre-defined optimal temperature values. Only works if length(zo) is S.
  # gam:       A vector of pre-defined mean dispersal distances. Only works if length(gam) is S. If no vector is specified, gamMean and gamSD are used to randomly generate the vector gam.
  # sig:       A vector of pre-defined thermal tolerance breadths. Only works if length(sig) is S. If no vector is specified, sigMean and sigSD are used to randomly generate the vector sig.
  # A:         Matrix of the relative competition coefficients between species. 
  # m:         A vector of pre-defined mortality probabilities (or a single value that will be shared for all species). These are probabilies, so the value must be between 0 and 1 (inclusive). Only works if length(m) is 1 or S.
  # gamMean:   Mean dispersal distance for randomized species. Default is based on Urban et al. 2012.
  # gamSD:     Standard deviation of dispersal distance for randomized species. Default is based on Urban et al. 2012.
  # sigMean:   Mean thermal tolerance breadth for randomized species. Default is based on Urban et al. 2012.
  # sigSD:     Standard deviation of thermal tolerance breadth for randomized species. Default is based on Urban et al. 2012.
  # lam:       Skewness in thermal tolerance. Default is based on Urban et al. 2012. (to have a mild decrease moving toward colder temperatures and a sharp decrease moving toward warmer temperatures).
  # B:         Area of integrated birth rate over all T for each species.
  # compType:  The type of competition in a string. Must be either "lottery" or "temp".
  
  # XW:        Window of analysis (to remove edge effects)
  
  # temp2d:    A matrix of pre-defined temperatures over x. Only works if nrow(temp2d) is L and ncol(temp2d) is W.
  # tempLow:   Lowest mean temperature on linear temperature gradient. temp1d(L)=tempLow
  # tempHigh:  Highest mean temperature on linear temperature gradient. temp1d(1)=tempHigh
  # tempRev:   If tempRev=T, then temp1d(1)=tempLow and temp1d(L)=tempHigh.
  # tempGH:    Hurst exponent for global temperature heterogeneity. H=0.5 is Brownian motion; 0.5<H<=1 is long-term positive autocorrelation and 0<=H<0.5 is long-term negative autocorrelation.
  # tempGSD:   Controls the magnitude of the global temperature heterogeneity.
  # tempLSD:   Standard deviation in temperature between microhabitats in each patch.
  # years:     Maximum number of years for initialization + climate change.
  # tempY:     A vector of pre-defined temperatures over time. Only works if length(tempY) is years.
  # tau:       Average temperature change per year.
  # tempYAC:   Temperature autocorrelation over time. Default is global temperature AC from 1880-1979.
  # tempYSD:   Temperature standard deviation over time. Default is global temperature SD from 1880-1979.
  # tempYXGH:  Hurst exponent for global heterogeneity in temperature change over time
  # tempYXGSD: Magnitude of global heterogeneity in temperature change over time
  # tempYXLSD: Standard deviation in temperature change between microhabitats in each patch
  
  # Q:         A matrix of pre-defined habitat quality over x. Only works if nrow(Q) is L and ncol(Q) is W.
  # QMean:     Average habitat quality for any microhabitat
  # QGH:       Hurst exponent for global habitat quality heterogeneity. H=0.5 is Brownian motion; 0.5<H<=1 is long-term positive autocorrelation and 0<=H<0.5 is long-term negative autocorrelation.
  # QGSD:      Controls the magnitude of the global habitat quality heterogeneity.
  # QLSD:      Standard deviation in habitat quality between microhabitats in each patch.
  
  
  
  
  ##########################################################
  # First, specify biological parameter values for all of the S species.
  
  # Optimal temperature for each species. These can be randomly picked from a uniform distribution or pre-defined.
  if(is.null(zo)){
    zo <- runif(S,9.9,30.1)
  } else {
    if(length(zo)!=S){
      stop("zo does not match the number of species!")
    }
  }
  
  # Dispersal distance for each species. These can be randomly picked from a lognormal distribution or pre-defined.
  if(is.null(gam)){
    # To use the lognormal distribution, we need to convert mean and SD values
    gamMu <- log(gamMean/sqrt(1+gamSD^2/gamMean^2))
    gamSig <- sqrt(log(1+gamSD^2/gamMean^2))
    gam <- rlnorm(S,gamMu,gamSig)
  } else {
    if(length(gam)!=S){
      stop("gam does not match the number of species!")
    }
  }
  
  # Thermal tolerance breadth for each species. These can be randomly picked from a lognormal distribution or pre-defined.
  if(is.null(sig)){
    # To use the lognormal distribution, we need to convert mean and SD values
    sigMu <- log(sigMean/sqrt(1+sigSD^2/sigMean^2))
    sigSig <- sqrt(log(1+sigSD^2/sigMean^2))
    sig <- rlnorm(S,sigMu,sigSig)
  } else {
    if(length(sig)!=S){
      stop("sig does not match the number of species!")
    }
  }
  
  # Competition coefficients between each pair of species species. By default, all coefficients are 1 (lottery competition), but A can be pre-defined.
  if(is.null(A)){
    A <- matrix(1,S,S)
  }else{
    if(!(nrow(A)==S & ncol(A)==S)){
      stop("A does not match the number of species!")
    }
  }
  
  # Yearly mortality probability for each species.
  if(any(m<0) | any(m>1)){
    stop("m should be between 0 and 1 (inclusive)")
  }
  if(length(m)==1){
    # If only one value is provided, all species will have the same mortality
    m <- rep(m,S)
  } else if(length(m)!=S){
    stop("m does not match the number of species!")
  }
  # Mortality is more convenient in this form
  M <- rep(m,W*L)
  
  # Make sure that compType is used correctly
  if(!(compType=="lottery" | compType=="temp")){
    stop("compType must be either 'lottery' or 'temp'!")
  }
  
  ##########################################################
  # Using the randomized species parameters, we derive other variables needed for computation.
  
  # zo helps define where the species' thermal optimum is, but mathematically this is not completely correct.
  # If we let zo be the true optimum, z is the value we plug into the reproduction function so that argmax(b)==zo 
  # We apply the function zAdjust to all values of zo to calculate z
  z <- mapply(zAdjust, sig, zo, lam, 2^13)
  
  # To speed up computation time, we define full dispersal kernels now.
  # The dispersal kernel uses q, a transformation of gam
  q <- sapply(1:S, function(i) 1+1/gam[i]-sqrt(1+1/gam[i]^2))
  l <- (-L):(L)
  k <- t(sapply(1:S, function(i) doubGeom(l,q[i])))
  k[k<10^(-15)] <- 0.0
  # K is a list of separate L by 2L dispersal kernel matrices
  # K[[s]] is the dispersal kernel matrix of species s
  # K[[s]][i,] is a vector of probabilities for a propagule in patch i to spread to patch j-L/2 (this is extra long to account for a propagule spreading beyond the limits of the ecosystem)
  K <- rep(list(matrix(0,L,2*L)),S)
  for(i in 1:S){
    Ki<-matrix(0,2*L,4*L)
    for(j in 1:(2*L)){
      Ki[j,j:(j+2*L)]<-k[i,]
      Ki[j,3/2*L]<-sum(Ki[j,1:(3/2*L)])
      Ki[j,5/2*L+1]<-sum(Ki[j,(5/2*L+1):(4*L)])
    }
    K[[i]]<-Ki[(L/2+1):(3*L/2),(3/2*L):(5/2*L)]
  }
  # Tolerance and reproductive strength have a tradeoff
  # Birth rate is adjusted so each species has roughly equal birth rate when integrated over all T
  # ro is a constant that adjusts to this reproductive output
  ro <- sapply(sig,function(x) rAdjust(x,B,lam,1e-06,2^13))
  
  ##########################################################
  # Put the biological parameters together into a single list, P
  P=list(S=S,
         z=z,
         gam=gam,
         sig=sig,
         lam=lam,
         A=A,
         M=M,
         ro=ro,
         zo=zo,
         K=K,
         compType=compType)
  
  ##########################################################
  # Next, we define environmental parameters.
  # Discrete spatial domain: from 1 to L (integers)
  x <- seq(1,L)
  # temp2d is the current temperature over all x and microhabitats
  if(is.null(temp2d)){
    
    temp1dr <- seq(tempHigh,tempLow,length=L)
    if(tempRev){
      temp1dr<-rev(temp1dr)
    }
    tempG <- tempVarH(L,tempGH)
    temp1d <- temp1dr+tempG*tempGSD/sd(tempG)
    if(W>1){
      temp2dr <- matrix(temp1d,L,W)
      tempLH <- matrix(rnorm(L*W,0,1),L,W)
      tempLH2 <- t(sapply(1:L, function(i) mean(tempLH[i,])+(tempLH[i,]-mean(tempLH[i,]))*tempLSD/sd(tempLH[i,])))
      tempLH3 <- t(matrix(sapply(1:L, function(x) sort(tempLH2[x,]-mean(tempLH2[x,]))),W,L))
      temp2d <- temp2dr+tempLH3
    } else{
      temp2d<-temp1d
    }
    
  } else {
    if(!(nrow(temp2d)==L & ncol(tempsd)==W)){
      stop("temp2d does not match environment size!")
    }
    temp1d<-rowMeans(temp2d)
  } 
  # tempY is a vector of the temperature over time
  if(is.null(tempY)){
    tempY<-0:years
    for (i in 1:years){
      # A new epsi is calculated for each time step
      tempY[i+1] <- tempYAC*tempY[i]+rnorm(1,0,tempYSD)*sqrt(1-tempYAC^2)
    }
  } else{
    if(length(tempY)!=years+1){
      stop("tempY does not match years!")
    }
  }
  # Tau is the temperature change over time in each patch and subpatch
  if(is.null(Tau)){
    tauLW <- matrix(tau*100,L,W)
    tauG <- tempVarH(L,tempYXGH)
    tauLW2 <- tauLW + tauG*tempYXGSD/sd(tauG)
    if(W>1){
      tauLH <- matrix(rnorm(L*W,0,tempYXLSD),L,W)
      tauLH <- t(matrix(sapply(1:L, function(x) tauLH[x,]-mean(tauLH[x,])),W,L))
    } else{
      tauLH <- 0
    }
    Tau <- tauLW2+tauLH
    Tau <- Tau/100
  } else{
    if(!(nrow(Tau)==L & ncol(Tau)==W)){
      stop("Tau does not match environment size!")
    }
  }
  
  # Habitat quality could differ in space or with species, but we will keep it constant for now
  if(is.null(Q)){
    # Q is updated, but it should only work if W==1
    QG <- tempVarH(L,QGH)
    Q <- (QMean*2)/(1+exp(-QG*QGSD/sd(QG)))
    Q[Q<=0.001]<-0.001
    Q2d <- array(Q,c(L,W))
    Qrep<-array(rep(Q2d,each=S),c(S,L,W))
  } else {
    if(!(length(Q)==L)){
      stop("Q does not match environment size!")
    }
  }
  
  ##########################################################
  # Put the abiotic parameters together into a single list, X
  
  X=list(L=L,
         x=x,
         XW=XW,
         temp1d=temp1d,
         tau=tau,
         Tau=Tau,
         Q=Q,
         Qrep=Qrep,
         tempY=tempY,
         tempYAC=tempYAC,
         tempYSD=tempYSD,
         temp2d=temp2d,
         W=W)
  
  ##########################################################
  # Export it all as a list
  
  return(list(P=P,X=X))
}

tempVarH <-  function(L,H,cZero=T){
  # This function adds some heterogeneity to the temperature gradient with fractional Brownian noise
  # See Keitt (2000)
  # Spectral representation of neutral landscapes
  # Landscape Ecology 15
  
  # H:      Hurst exponent (should be between 0 and 1). It relates to the autocorrelation.
  ###        When H is near 1, this function has positive long-term positive autocorrelation and will look relatively smooth.
  ###        When H=0.5, this function is Brownian motion.
  ###        When H is near 0, then autocorrelation is negative and positive values will more often be followed by negative values (and vice versa).
  # L:     Length of the temperature gradient.
  # cZero: If T, it will center the whole output so the mean is 0.
  
  
  # random phases uniformly distributed on [0,2pi]
  phif <- runif(L)*2*pi
  
  # adjusted exponent for amplitudes
  betaH <- 1+2*H
  
  # uniformly distributed random numbers
  xf <- rnorm(L)
  # to form the amplitudes
  af <- 1/seq(1,L)^(betaH/2)*xf
  
  # complex coeffcients
  cf <- af*exp(1i*phif)
  
  #  real part of the inverse fourier transform
  tH <- Re(ifft(cf))
  
  # center it around zero?
  if(cZero){
    tH <- tH-mean(tH)
  }
  
  # multiply the output to increase the magnitude of the heterogeneity
  # add that to the the temperature gradient
  return(tH)
}


doubGeom<-function(x,q){
  # Probability mass function for "double geometric" distribution
  # x: distance from origin to landing spot
  # q: probability of remaining in a given patch (and not continuing to move); see supplemental
  
  return((q/(2-q)*(1-q)^abs(x)))
}

rAdjust<-function(sig,B,lam=-2.7,eps=1e-06,len=2^13){
  # This function creates a constant to adjust the reproduction rate so that the area under the curve is roughly equal for all species
  # sig: Thermal tolerance width of a species
  # B:   Desired total integrated area of positive growth
  # lam: Skewness in thermal tolerance
  # eps: Precision of estimate
  # len: Length of temperature vector. Higher values are more precise
  
  # Set up an extended version of a linear tempereature gradient
  temp <- seq(-100,100,length=len)
  
  # The actual optimal temperature is not important here, so we use the center of the temperature gradient
  z <- 20
  r <- exp(-(temp-z)^2/sig^2)*(1+erf(lam*(temp-z)/sig))-1
  
  bL <- -125; bH <- 125
  
  # Binary search for a value of ro such that exp(ro*r) integrates to B over all temperature values where exp(ro*r) is positive
  
  for(i in 1:500){
    bM <- (bL+bH)/2
    R <- exp(bM*r)
    G <- trapz(temp,(R-1)*(R>1))
    differ<- G-B
    if(abs(differ)<eps){
      break
    } else{
      if(differ>0){
        bH<-bM
      } else{
        bL<-bM
      }
    }
  }
  return(bM)
}

zAdjust<-function(sig,zo,lam=-2.7,L=2^13){
  # The reproduction function in Urban et al. 2012 is useful for creating the shape of the reproduction rate over temperature
  # However, the z_i "optimal temperature" doesn't end up where we might expect it to be
  # This function adjusts so that argmax_{temp1d}(R_i)=z_i
  # sig: The thermal tolerance width of a species
  # z:   Optimal temperature of species
  # lam: Skewness in thermal tolerance
  # len: Length of temperature vector. Higher values are more precise
  
  # Set up an extended version of a linear tempereature gradient
  temp <- seq(-100,100,length=L)
  
  # We need to calculate the difference between the expected optimal temperature and the actual optimal temperature
  # To do so, we begin with a baseline at zc=20
  zc <- 20
  
  # Calculate the baseline reproductive rate
  r<-exp(-(temp-zc)^2/sig^2)*(1+erf(lam*(temp-zc)/sig))-1
  
  # index for which temperature has the maximum reproductive output with the baseline
  iZ<-which.max(r)
  # index for baseline optimal temperature
  oZ<-which.min(abs(temp-zc))
  # index for desired optimal temperature
  tZ<-which.min(abs(temp-zo))
  
  # adjusted z to make optimal temperature in the right place
  z<-temp[tZ+oZ-iZ]
  
  return(z)
}

commSimulate <- function(n,P,X,y=1,years=100,init=F,extInit=F,extThresh=100,AM=NULL){
  # This simulates a community, n, over yars
  # n:         Initial population sizes. SxLxW array of population nonnegative integers.
  # P:         List of biotic variables
  # X:         List of abiotic variables
  # years:     How many time steps to run the model.
  # extInit:   If T, the simulation stops running after extThresh time steps without any extinctions. When attempting to initialize a stable community, consider setting extInit to T.
  # extThresh: If extInit==T, then the simulation stops once extInit time steps have passed without any extinctions.
  # manage:    A list with a bunch of management options. If left blank, the model generates a list where all management options are set to FALSE.
  
  # Make an all-FALSE management list if none is provided
  if(is.null(AM)){
    AM <- amSetup(P,X)
  }
  if(init==T){
    Tau <- 0
  } else{
    Tau <- X$Tau
  }
  
  # First, we set up a matrix to save the total population size over time
  N <- matrix(0,P$S,years+1)
  # Record the total initial population size of each species across the whole ecosystem
  N[,1] <- apply(n,1,sum)
  
  # Temperature changes over time, so we need to adjust this over the course of the model
  temp2d0 <- X$temp2d
  
  # For output, we want to keep track of the average temperature over time, and we do that with temps
  temps <- seq(0,years)
  temps[1] <- mean(temp2d0+X$tempY[y])
  
  # Keep track of tempY and tau outside of X
  tempY <- X$tempY
  
  # Run the model for a number of time steps equal to 'years'
  for (i in 1:years){
    # Temperature changes before each time step
    X$temp2d=temp2d0+Tau*i+tempY[i+y]
    X$temp1d=rowMeans(X$temp2d)
    # save the mean temperature
    temps[i+1]=mean(X$temp2d)

    if(sum(N[,i])>0){
      ts <- timeStep(n,P,X,N[,i],i,AM)
      n <- ts$n
      AM <- ts$AM
    }
    # Record the population size
    N[,i+1]<-apply(n,1,sum)
  }
  return(list(n=n,N=N,temps=temps,AM=AM))
}

timeStep <- function(n,P,X,N,t,AM){
  # Cycle through each step of the model.
  # Each time step could be one "year" or one "generation", but ultimately it runs through each part of the life cycle in an order determined by lcOrder.
  # The various management techniques can optionally be added to the model between any two of the required steps
  # reproduction -> dispersal -> density dependence
  # Reproduction
  n1 <- reproduce(n,X$L,X$W,P$S,P$z,P$sig,P$ro,P$lam,X$temp2d)

  
  # If assisted migration is occurring in this simulation
  if(AM$AM==T){
    am <- assistMigrate(n1,N,X$L,X$W,X$Q,X$temp1d,X$tau,X$XW,P$S,t,AM)
    # These are the individuals (of all species) that will still be dispersing normally
    n1 <- am$nD
    # These are the individuals that were relocated
    nR <- am$nR
    # Update the time since last relocation vector
    AM$tLR <- am$tLR
    # Make note of the times when the species was relocated
    AM$relTimes[am$SReloc,t] <- 1
    # Subtract the moves from movesLeft
    AM$movesLeft<-AM$movesLeft-length(am$SReloc)
  } else{
    nR <- n1*0
  }
  
  # To save computation time, we don't need to use the dispersal function on species without any propagules
  # dS are the species with extant propagules
  dS<-which(rowSums(n1)>0)
  # Thus, as long as there is at least one species dispersing, go through with dispersal on dS species
  if(!isempty(dS)){
    # Slice the array so we are only using those that will be dispersing
    dispn <- n1[dS,,,drop=F]
    # Now run the disperse function on the slice
    dispn2 <- disperse(dispn,X$L,X$W,P$S,P$K[dS])
    # Preallocate the dispersed array
    n2 <- n1*0
    # Add the dispersed individuals to that array
    n2[dS,,] <- dispn2
  } else {
    # If there are no propagules at all, n2 is just n1
    n2<-n1
  }
  
  n2 <- n2+nR
  # The survive function simulates mortality of adults
  n3 <- n2+survive(n,X$L,X$W,P$S,P$M)
  
  # All inidividuals then compete
  # (At this point, adults and offspring have equal competitive ability. We could change the compete function if this should change.)
  n4 <- compete(n3,X$L,X$W,P$S,X$Qrep,P$A,P$compType,P$z,P$sig,P$ro,P$lam,X$temp2d)
  ts<-list(n=n4,AM=AM)
  return(ts)
}

bi <- function(z,sig,ro,lam,temp){
  # reproductive rate function
  op<-ro*(exp(-((temp-z)/sig)^2)*(1+erf(lam*(temp-z)/sig))-1)
  return(op)
}

reproduce <- function(n,L,W,S,z,sig,ro,lam,temp2d){
  # The number of offspring born for each species in each location is a Poisson random variable with mean r*n
  
  # The base reproductive rate is a skewed function, adjust such that min(r)=0 and max(r)=2
  # Each species will have a different reproductive rate depending on the temperature at that space.
  # r is the "continuous" form of the birth rate
  r <- sapply(1:S, function(i) bi(z[i],sig[i],ro[i],lam,temp2d))
  # R turns it into a discrete form
  R <- exp(r)
  # This just turns it into the correct array format
  R <- aperm(array(R,c(L,W,S)),c(3,1,2))
  
  # Mean number of offspring
  rn <-c(R*n)
  
  # The number of offspring is a Poisson random variable with mean=R*n
  nr<-array(sapply(rn, function(x) rpois(1,x)),c(S,L,W))
  
  return(nr)
}

disperse <- function(n,L,W,S,K){
  # Each individual spreads throughout the spatial landscape with a random double geometric dispersal kernel determined by the species' mean dispersal distance, gam[i].
  # For each species in each location, disperse!
  
  # Si is the total number of species that are dispersing 
  Si <- nrow(n)
  
  if(is.null(Si)){
    # When there is 1 species, this function gets confused, so we can fix it here
    n <- array(n,c(S,L,W))
    Si <- S
  }
  
  # Flatten the metapopulation so that all microhabitats in one patch are summed together
  # This makes n1 an SxL matrix
  n1 <- apply(n,c(1,2),sum)
  # This disperses the propagules with multinomial random vectors across the temperature gradient X
  n2 <- t(sapply(1:Si, function(j) disperseMulti(n1[j,],L,K[[j]])))
  # This distributes the propagules randomly into the microhabitats for each patch x
  n3 <- c(sapply(1:Si, function(i) t(sapply(1:L,function(j) rebin(sample(1:W,n2[i,j],replace=T),W)))))
  # This just reformats n3 so it is in the previous SxLxW form
  n4 <- aperm(array(n3,c(L,W,Si)),c(3,1,2))
  
  # And now it's ready to go
  return(n4)
}

disperseMulti <- function(n,L,K){
  # Used in in the disperse function
  # To save computation time, we only use the mulinomial random number generator for patches where local n is positive
  # y is just there to mark which indices we are going to run the multinomial generator
  y <- which(n>0)
  # Run the multinomial random generator
  n1 <- sapply(y,function(x) rmultinom(1,n[x],K[x,]))
  
  # Now we add all of these vectors together
  if(length(y)>1){
    # (Assuming that propagules dispersed from more than one patch)
    n2 <- rowSums(n1)
  } else{
    # (Otherwise, no summation is necessary)
    n2 <- n1
  }
  # This just cuts off the edges that are removed from the model (since we have absorbing boundaries)
  n3 <- n2[2:(L+1)]
  return(n3)
}

survive <- function(n,L,W,S,M){
  # Adults survival is stochastic
  # First, we flatten n so it is one long vector (SLWx1) (just like M)
  nv <- c(n)
  
  # We can save computation time if we skip over cases where all species have 0 or 1 mortality probabilities.
  if(all(M==1)){
    # When all M is 1, all adults die
    ns <- n*0
  } else if(all(M==0)){
    # When all M is 0, all adults live
    ns <- n
  } else{
    # To save computation time, we find out which patches have living adults with some probability of surviving
    # wMort are all of the patches that fit this
    wMort <- which(nv>0 & M<1)
    
    # nsv is the full vector form of surviving adults. Pre-allocated to 0 for all.
    nsv <- 0*nv
    
    # Assuming that at least one patch with adults that might survive, calculate stochastic survival
    if(!isempty(wMort)){
      nsi<-sapply(wMort, function(i) rbinom(1,nv[i],1-M[i]))
      nsv[wMort]<-nsi
    }
    # Convert the output into original SxLxW form
    ns<-array(nsv,c(S,L,W))
  }
  
  # ns is all of the surviving adults
  return(ns)
}


compete <- function(n,L,W,S,Qrep,A,compType='lottery',z=NULL,sig=NULL,ro=NULL,lam=NULL,temp2d=NULL){
  # The density dependence in this model is roughly a Beverton-Holt model that includes both interspecific and intraspecific competition
  # Each individual has a random chance of survival based on a variety of conditions
  
  # Competition coefficients depend on interactions between each species and the temperature at the location at the time
  # These can be thought of as temperature-varying Lotka-Volterra competition coefficients
  # Probability of survival depends on competition coefficients, number of individuals of each different species at that location, and the quality of the habitat at that location
  # Competition works differently dependenting on whether it is temperature-dependent or pure lottery competition
  if(compType=="temp"){
    # Use the same reproduction temperature dependence to determine competitive pressure
    r <- sapply(1:S, function(i) bi(z[i],sig[i],ro[i],lam,temp2d))
    R <- aperm(array(exp(r),c(L,W,S)),c(3,1,2))
  } else if(compType=="lottery"){
    # All individuals are equal
    R <- 1
  }
  
  # Convert Q into an SxLxW array
  # QR determines habitat quality the species
  QR <- 1/(Qrep*R)
  # nR is used to determine species interactions
  nR <- R*n
  # This puts the species interactions together
  anR <- sapply(1:S, function(s) colSums(A[s,]*nR))
  anR <- aperm(array(anR,c(L,W,S)),c(3,1,2))
  
  # Convert this into survival probability
  p <- 1/(1+QR*anR)
  
  # Binomial random variables to see who survives
  nc <- (sapply(1:S,function(s) mapply(rbinom,1,c(n[s,,]),p[s,,])))
  # Converted into proper SxLxW form
  nc2 <- array(t(nc),c(S,L,W))
  
  return(nc2)
}

assistMigrate<-function(n,N,L,W,Q,temp1d,tau,XW,S,t,AM){
  # Attach AM parameters
  tLR <- AM$tLR; targs <- AM$targs; eta <- AM$eta; tCD <- AM$tCD; movesLeft <- AM$movesLeft
  
  # Preallocate the output arrays
  # nR is the array of relocated individuals
  nR <- n*0
  # nD is the array of individuals that will disperse naturally instead of assisted migration
  nD <- n
  
  ##########################################################
  # DO WE DO ANY ASSISTED MIGRATION DURING THIS TIME STEP? #
  ##########################################################
  # Which species need to be relocated?
  # Must be a target species, not during a cooldown period, and population less than the threshold but greater than 0
  SReloc <- which((1:S)%in%targs & (t-tLR)>tCD & N<eta & N>0 & movesLeft>0 & rowSums(flatComm(n)[,XW])>0)
  SRL <- length(SReloc)
  
  # If there are more species to move than there are remaining moves, randomly select species to move  
  if(SRL>movesLeft){
    SReloc <- sample(SReloc,movesLeft)
    SRL <- movesLeft
  }
  
  # We only need to go through this bit if there are going to be any relocations during this time step
  if(!isempty(SReloc)){
    
    # Attach AM parameters
    rho <- AM$rho; mu <- AM$mu; zEst <- AM$zEst; xLoc <- AM$xLoc; recRad <- AM$recRad; lowQualThresh <- AM$lowQualThresh
    # Preallocate the relocation array
    nXWReloc <- n[SReloc,,,drop=FALSE]*0
    
    # Because relocation is occuring, we can update the tLR (time of last relocation) vector
    tLR[SReloc] <- t
    
    # We don't need to worry about microhabitats here, so we can flatten out the population array a bit
    # This makes n1 an SxL matrix
    # To help with this, set up the populations for each SReloc species into a vector for all patches with microhabitats summed up (SReloc x L)
    nv <-  flatComm(n[SReloc,,,drop=F])
    
    # Population size of selected species
    NSProps <- rowSums(nv)
    ##################################################
    # HOW MANY AND WHICH INDIVIDUALS DO WE RELOCATE? #
    ##################################################
    # NDon is the total number of donated propagules for each species in SReloc
    NDon <- ceil(NSProps*rho)
    
    # Now we pick the actual individuals out of the metapopulations
    # We can pick them in different ways
    # 2 is picking individuals from the trailing edge
    # First, we preallocate the nDon
    nDon <- nv*0
    for(s in 1:SRL){
      # We convert the spatial distribution of local populations into a vector that just shows where each individual is
      nUnbin <- unbin(nv[s,])
      # Pick the first NDon[,s] on the trailing edge
      nUnbinDon <- nUnbin[1:NDon[s]]
      # and put it back into the regular format
      nDon[s,] <- rebin(nUnbinDon,length(nv[s,]))
    }
    
    # Remove the donated indivudals from the nv array
    nvDisp <- nv-nDon
    # Preallocate an array for relocated individuals
    nvReloc <- nv*0
    
    ############################
    # WHO SURVIVES RELOCATION? #
    ############################
    # Total individuals surviving
    NDonS <- sapply(1:SRL, function(s) rbinom(1,NDon[s],mu))
    
    #########################
    # WHERE DO WE RELOCATE? #
    #########################
    # Which patch is closest to the estimated thermal optimum + xLoc patches ahead
    locOpt <- rep(NA,SRL)
    for(i in 1:SRL){
      s<-SReloc[i]
      loc <- which( (temp1d+tau*xLoc)<zEst[s] & Q>=lowQualThresh)
      if(length(loc)==0){
        locOpt[i] <- length(temp1d)
      } else{
        locOpt[i] <- min(loc)
      }
    }
    # If locOpt is too big or too small, we need to fix that
    for(i in 1:SRL){
      if(locOpt[i]<(1+recRad)){
        locOpt[i] <- 1+recRad
      }else if(locOpt[i]>(L-recRad))
        locOpt[i] <- L-recRad
    }
    
    # Identify the locations that receive inidividuals
    locOptAll <- t(sapply(1:SRL, function(s) (locOpt[s]-recRad):(locOpt[s]+recRad)))
    
    # Relocate the individuals based on shape
    # 1 is a square shape
    for(i in 1:SRL){
      # Length of recipient location
      recSize<-(2*recRad+1)
      # Evenly distribute individuals through the recipient location
      nvReloc[i,locOptAll[i,]] <- nvReloc[i,locOptAll[i,]]+floor(NDonS[i]/recSize)
      # There will probably be a few extras after all is even. Place these randomly.
      extras <- sample(locOptAll[i,],NDonS[i]%%recSize)
      nvReloc[i,extras] <- nvReloc[i,extras]+1
      
      # Now we redistribute these amongst the width of the ecosystem
      # Similarly to above
      for(x in locOptAll[i,]){
        nXWReloc[i,x,] <- floor(nvReloc[i,x]/W)
        extras <- sample(1:W,nvReloc[i,x]%%W)
        nXWReloc[i,x,extras] <- nXWReloc[i,x,extras]+1
      }
    }
    # Now put the relocated and dispersing individuals into full arrays
    nR[SReloc,,] <- nXWReloc
    # The non-relocated individuals will disperse naturally
    # When running the disperse() function, all microhabitats are combined, so there's no reason to separate these into microhabitats
    nD[SReloc,,1] <- t(nvDisp)
  }
  # Ultimately, we want to output the relocated array, the dispersing array, the species that were relocated, and the updated time since last relocation
  output <- list(nR=nR,
                 nD=nD,
                 SReloc=SReloc,
                 tLR=tLR)
  return(output)
}

unbin <- function(v){
  # Convert vector of population sizes over x into a vector of the location for each individual
  L<-sum(v)
  ub<-matrix(0,L)
  j<-0
  for(i in which(v>0)){
    ub[(j+1):(j+v[i])]<-i
    j<-j+v[i]
  }
  return(c(ub))
}

rebin <- function(ub,L){
  # Converts a vector of individual locations into a vector of population sizes over x
  v<-1:L
  for(i in 1:L){
    v[i]<-length(which(ub==i))
  }
  return(v)
}


divInd <- function(n,type="alpha",index="invSimp",aMean=T){
  # Calculates the diversity index of a matrix of population sizes
  if(type=='alpha'){
    p <- t(n)/colSums(n)
    p2 <- rowSums(p^2)
    p2[is.nan(p2)] <- 1
    
    if(aMean==T){
      if(index=="invSimp"){
        D <- mean(1/p2)
      } else if(index=="giniSimp"){
        D <- 1- mean(p2)
      } else if(index=="richness"){
        D <- mean(colSums(n>0))
      }
    } else{
      if(index=="invSimp"){
        D <- 1/p2
      } else if(index=="giniSimp"){
        D <- 1- p2
      } else if(index=="richness"){
        D <- colSums(n>0)
      }
    }
  } else if(type=='gamma'){
    p <- rowSums(n)/sum(n)
    if(index=="invSimp"){
      D <- 1/sum(p^2)
    } else if(index=="giniSimp"){
      D <- 1- sum(p^2)
    } else if(index=="richness"){
      D <- sum(rowSums(n>0)>0)
    }
  }
  return(D)
}

commTrim <- function(n,P,X){
  # Remove extinct species from n and P
  nFlat <- t(sapply(1:P$S, function(s) rowSums(matrix(n[s,,,drop=FALSE],X$L,X$W))))
  extant <- which(rowSums(nFlat)>0)
  P$S <- length(extant)
  P$z <- P$z[extant,drop=FALSE]
  P$gam <- P$gam[extant,drop=FALSE]
  P$sig <- P$sig[extant,drop=FALSE]
  P$A <- P$A[extant,extant,drop=FALSE]
  P$ro <- P$ro[extant,drop=FALSE]
  P$zo <- P$zo[extant,drop=FALSE]
  P$K <- P$K[extant]
  n <- n[extant,,,drop=FALSE]
  ct <- list(n=n,P=P)
  return(ct)
}

flatComm <- function(n){
  # Take a community with subpopulations and flatten it into a single patch each
  nf <- apply(n,c(1,2),sum)
  return(nf)
}

qRange <- function(nFi,q1=0.025,q2=0.975,type=7){
  # Find the range of a species from quantiles (to cut of long tails that aren't representative of true range)
  # To find full range, you can set q1=0 and q2=1
  # nFi is the flattened (1d) vector of population size over space
  # q1 is the low quantile
  # q2 is the high quantile
  # type determines how to find the quantile. 7 is the default.
  if(sum(nFi)<=0){
    r<-0 # if the total population size is 0 (or somehow less?) the range is 0
  }else{ # otherwise, it is determined from quantiles over space
    ubn <- unbin(nFi) # convert the population matrix to individual location over space
    q <- quantile(ubn,c(q1,q2),names=F,type=type) # find the high and low quantiles 
    r <- q[2]-q[1] # find the difference. this is the range
  }
  return(r)
}


# Restoration methods


amSetup <- function(P,X,years=100,
                    AM=NULL,
                    movesLeft=50,
                    eta=50,tCD=5,rho=0.55,mu=0.8,zEst=P$zo,recRad=2,lowQualThresh=5,xLoc=5,
                    lowTempCutoff=14.86){
  if(is.null(AM)){
    # If you have no plans for AM at all, there is no need to define anything but AM$AM.
    AM <- list(AM=F)
  } else{
    targs <- which(P$zo>lowTempCutoff)
    if(length(targs)==0){
      AM<-list(AM=F)
    } else{
      tLR<-rep(0,P$S)
      relTimes<-matrix(0,P$S,years)
      AM <- list(AM=AM, movesLeft=movesLeft, targs=targs, eta=eta, tCD=tCD, rho=rho, mu=mu, zEst=zEst, recRad=recRad, lowQualThresh=lowQualThresh, xLoc=xLoc, tLR=tLR, relTimes=relTimes)
    }
  }
  return(AM)
}


corridors <- function(Q,effort,eps=0.0001,trys=1000){
  # This runs a searching algorithm to find a QNew that increases the minimum value of Q throughout the entire landscape so that sum(QNew)==sum(Q)+effL
  # Q is the habitat quality (carrying capacity) over space
  # effort is the amount of increase in carrying capacity over all space in units of L (length of space). Raw increase will be effort*L
  # eps is the amount of precision in finding this new Q (smaller is more precise)
  # trys is the maximum number of times we loop through the search
  L<-length(Q) # length of space
  effL<-effort*L # raw effort value
  QSum<-sum(Q) # total of all Q from the old vector
  QSumEff<-QSum+effL # sum of all Q in target vector
  QLL<-min(Q)  # high value for the new minimum of Q
  QLH<-max(Q)+effort  # low value for the new minimum of Q
  # narrow down during each step in the loop
  for(i in 1:trys){ 
    QLM <- (QLL+QLH)/2 # middle vale for the new minimum of Q
    QM <- Q # form a new candidate Q
    QM[Q<QLM]<-QLM # everything below QLM is now WLM
    QSumM<-sum(QM) # what's the sum of this new Q?
    QDiffM<-QSumM-(QSumEff) # difference between old Q and new candidate Q + effL
    if(abs(QDiffM)<eps){
      break # if the absolute difference is less than eps, we've found the new Q
    } else{ # if not, we update high and low values for the new minimum
      if(QDiffM>0){
        QLH<-QLM # if the sum of QM is higher than the sum of Q + effL, it becomes the new maximum
      } else{
        QLL<-QLM # otherwise it is the new minimum
      }
    }
  }
  QNew <- QM # QNew is the corridor-restored value of Q
  return(QNew)
}

steppingStone <- function(Q,effort,islandSize=1/2,eps=0.0001,trys=1000){
  Qc<-corridors(Q,effort,eps,trys)
  Qi<-Q
  Qdiff<-Qc-Q
  QLocs<-which(Qdiff>0)
  group <- c(0, cumsum(diff(QLocs) != 1))
  groupInfo <- tapply(QLocs, group, summary)
  for(i in 1:length(groupInfo)){
    xH<-groupInfo[[i]][6]
    xL<-groupInfo[[i]][1]
    area<-sum(Qdiff[xL:xH])
    island<-islandLocs(xL,xH,islandSize)
    Qi[island]<-Qi[island]+area/length(island)        
  }
  QNew<-Qi
  return(QNew)
}

islandLocs <- function(gL,gH,islandSize){
  gapSize <- gH-gL
  gM <- round(gL+gapSize/2)
  iL <- round(gM-gapSize*islandSize/2)
  iH <- iL+islandSize*gapSize-1
  island <- iL:iH
  return(island)
}

allRestore <- function(Q,effort){
  QNew <- Q+effort
  return(QNew)
}

reinforce <- function(Q,effort,cutoff=0.75){
  L <- length(Q)
  #find out the locations with habitat quality higher than cutoff quantile
  QLocs<- which(Q>=quantile(Q,cutoff))
  #calculate the restoration effort for each individual high-quality patch 
  effortX <- effort*L/length(QLocs)
  
  QNew <- Q
  QNew[QLocs] <- Q[QLocs]+effortX
  
  return(QNew)
}

expand <- function(Q,effort,cutoff=0.75){
  thresh <- quantile(Q,cutoff)
  Qs1 <- Q
  Qs2 <- Q
  QSumTarg <- sum(Q+effort)
  
  for(i in 1:length(Q)){
    QLocsL <- which(Qs2<thresh)
    if(isempty(QLocsL)){
      QNew <- corridors(Q,effort)
      break
    }
    group <- c(0, cumsum(diff(QLocsL) != 1))
    groupInfo <- tapply(QLocsL, group, summary)
    lengths <- table(group)
    for(j in 1:length(groupInfo)){
      xH<-groupInfo[[j]][6]
      xL<-groupInfo[[j]][1]
      Qs2[c(xL,xH)]<-thresh
    }
    if(sum(Qs2)>QSumTarg){
      QsDiff <- Qs2-Qs1
      mult <- (QSumTarg-sum(Qs1))/sum(QsDiff)
      QNew <- Qs1+mult*QsDiff
      break
    } else{
      Qs1 <- Qs2
    }
  }
  return(QNew)
}

changeHeterogeneity <- function(temp2d,SD2){
  SD1 <- mean(sapply(1:512, function(x) sd(temp2d[1,])))
  temp2d2 <- rowMeans(temp2d)+(temp2d-rowMeans(temp2d))*(SD2/SD1)
  return(temp2d2)
}

QChange <- function(X,Q,S){
  X$Q<-Q
  Q2d <- array(Q,c(X$L,X$W))/X$W
  X$Qrep<-array(rep(Q2d,each=S),c(S,X$L,X$W))
  return(X)
}


topSpecies <- function(n){
  L<-nrow(n[1,,]) # Length (number of patches)
  W<-ncol(n[1,,]) # Width (number of subpatches per patch)
  S<-nrow(n[,,1]) # Number of species
  tSpec <- matrix(NA,L,W) # Output matrix
  for(i in 1:L){
    for(j in 1:W){
      nij <- n[,i,j] # All species' populations in the subpatch
      if(sum(nij)>0){
        nij<-nij+runif(S,0,0.001) # Add a small amount of randomness to break ties
        tSpec[i,j]<-which.max(nij) # Which species has the highest population size
      }
    }
  }
  return(tSpec)
}

vComTop<- function(n,S=NULL){
  if(is.null(S)){
    S<-length(unique(c(tSpec[!is.na(tSpec)]))) # We only need colors for species that are represented
  }
  tSpec<-topSpecies(n) # Create top view matrix
  image(tSpec,col=rainbow(S))
}

vComSide<-function(n,lwd=2,ymax=max(n),S=NULL){
  if(is.null(S)){
    S<-nrow(n[,,1]) # All species with a row are represented
  }
  nFlat<-sapply(1:S, function(s) rowSums(n[s,,])) # Flatten the subpatches in the community
  matplot(nFlat,type="l",lwd=lwd,lty=1,ylim=c(0,ymax),col=rainbow(S))
}

vComTime <- function(N,year1=2,year2=NULL,yearAd=0){
  S<-nrow(N) # Number of species
  if(is.null(year2)){
    year2<-ncol(N) # If year2 is not specified, it is determined by the size of the matrix
  }
  matplot((1:(year2-year1+1)),t(N[,year1:year2]),type="l",lwd=2,lty=1,ylim=c(0,max(N)),col=rainbow(S))
}

commOut <- function(id,tempLSD,tempYSD,sim1,sim2,X,P,Qtype,treat,effort,M1,M2){
   # relevant community data
  comm <- id
  n1 <- flatComm(sim1$n); n2<-flatComm(sim2$n)
  temps1 <- sim1$temps; temps2 <- sim2$temps
  tau <- X$tau; gam <- P$gam; sig <- P$sig; zo <- P$zo
  if(treat>6){
    relocs <- sum(sum(sim2$AM$relTimes))
  } else{
    relocs <- 0
  }
  S<-P$S
  n1M <- n1[,M1]; n2M <- n2[,M2]
  surv1 <- which(rowSums(n1M)>0); surv2 <- which(rowSums(n2M)>0)
  S1 <- sum(rowSums(n1)>0); S2 <- sum(rowSums(n2)>0)
  S1M <- sum(rowSums(n1M)>0); S2M <- sum(rowSums(n2M)>0)
  div1 <- c(divInd(n1M,"alpha","invSimp"),divInd(n1M,"alpha","richness"),divInd(n1M,"gamma","invSimp"),divInd(n1M,"gamma","richness"))
  div2 <- c(divInd(n2M,"alpha","invSimp"),divInd(n2M,"alpha","richness"),divInd(n2M,"gamma","invSimp"),divInd(n2M,"gamma","richness"))
  totN1 <- sum(n1M); totN2 <- sum(n2M)
  range1 <- sapply(surv1, function(s) qRange(n1[s,])); range2 <- sapply(surv2, function(s) qRange(n2[s,]))
  r1 <- mean(range1); r2 <- mean(range2)
  disp1 <- c(mean(gam[surv1]),sd(gam[surv1])); disp2 <- c(mean(gam[surv2]),sd(gam[surv2]))
  toler1 <- c(mean(sig[surv1]),sd(sig[surv1])); toler2 <- c(mean(sig[surv2]),sd(sig[surv2]))
  opt1 <- c(mean(zo[surv1]),sd(zo[surv1])); opt2 <- c(mean(zo[surv2]),sd(zo[surv2]))
  ext2 <- sum(rowSums(n2[surv1,,drop=F])==0)
  temp1 <- c(temps1[length(temps1)]-temps1[1], sd(temps1)); temp2 <- c(temps2[length(temps2)]-temps2[1],sd(temps2-seq(20,20+tau*(length(temps2)-1),length=length(temps2))))
  output <- c(comm,tempLSD,tempYSD,Qtype,treat,effort,relocs,S,
              S1,S2,div1,div2,totN1,totN2,
              r1,r2,disp1,disp2,toler1,toler2,opt1,opt2,temp1,temp2,ext2)
  return(output)
}

specOut <- function(id,s,tempLSD,tempYSD,sim1,sim2,X,P,Qtype,treat,effort,M1,M2){
  comm <- commOut(id,tempLSD,tempYSD,sim1,sim2,X,P,Qtype,treat,effort,M1,M2)
  
  n1 <- flatComm(sim1$n); n2<-flatComm(sim2$n)
  n1M <- n1[,M1]; n2M <- n2[,M2]
  temps1 <- sim1$temps; temps2 <- sim2$temps
  tau <- X$tau; gam <- P$gam; sig <- P$sig; zo <- P$zo
  if(treat>6){
    relocs <- sum(sum(sim2$AM$relTimes))
  } else{
    relocs <- 0
  }
  S<-P$S
  N1<-sim1$N; N2<-sim2$N
  t1<-ncol(N1); t2<-ncol(N2)
  alive1 <- which(rowSums(n1)>0)
  status1 <- (sum(n1[s,])>0)+(sum(n1M[s,])>0); status2 <- (sum(n2[s,])>0)+(sum(n2M[s,])>0)
  ext2 <- (status1==2 & status2==0)*1
  if(ext2){
    ttExt2 <- min(which(N2[s,]==0))-1
  } else{
    ttExt2 <- NA
  }
  
  zF <- zo[s]
  dispF <- gam[s]
  tolerF <- sig[s]
  NFi <- N1[s,t1]; NFf <- N2[s,t2]; NFMin <- min(N2[s,]); NFMax <- max(N2[s,]); NFMean <- mean(N2[s,]); NFSD <- sd(N2[s,])
  
  rF1 <- qRange(n1[s,]); rF2 <- qRange(n2[s,])
  
  zA1 <- zo[alive1]

  sPole1 <- which(zA1<zF); sEquator1 <- which(zA1>zF)
  if(isempty(sPole1)){
    zP1 <- NA; dispP1 <- NA; tolerP1 <- NA; NP1 <- NA; rP1 <- NA 
  } else{
    P<-alive1[sPole1[which.min(zF-zA1[sPole1])]]
    zP1 <- zF-zo[P]; dispP1 <- gam[P]; tolerP1 <- sig[P]; NP1 <- N1[P,t1]; rP1 <- qRange(n1[P,])
  }
  if(isempty(sEquator1)){
    zE1 <- NA; dispE1 <- NA; tolerE1 <- NA; NE1 <- NA; rE1 <- NA 
  } else{
    E<-alive1[sEquator1[which.min(zA1[sEquator1]-zF)]]
    zE1 <- zo[E]-zF; dispE1 <- gam[E]; tolerE1 <- sig[E]; NE1 <- N1[E,t1]; rE1 <- qRange(n1[E,])
  }
  output <- c(comm,s,status1,status2,ext2,ttExt2,
              NFi,NFf,NFMin,NFMax,NFMean,NFSD,
              rF1,rF2,zF,dispF,tolerF,
              NP1,rP1,zP1,dispP1,tolerP1,NE1,rE1,zE1,dispE1,tolerE1)
  return(output)
}


simulation <- function(id){
  set.seed(id)
  
  tempYSD <- runif(1,0,1)
  tempLSD <- runif(1,0,2)
  
  S <- 64     # number of species
  L <- 512     # length in space
  W <- 8       # width of space

  iYears<-500  # number of years to initialize the model  Q 
  ccYears<-100 # number of years to run climate change after initialization
  
  
  # Set up standard quality vectors
  #Y1<-rep(66,512)
  #Y2<-round(c(rep(66,64),64*(1+sin((1:400)*pi/100))+2,rep(66,48)))
  Y3<-round(c(rep(66,64),64*(1+sin((1:400)*pi/50))+2,rep(66,48)))
  #Y4<-round(c(rep(66,64),64*(1+sin((1:400)*pi/25))+2,rep(66,48)))
  Y5<-round(c(rep(66,64),64*(1+sin((1:400)*pi/12.5))+2,rep(66,48)))
  
  # Put them all in a matrix
  Q<-matrix(c(Y3,Y5),L,2)
  
  M1 <- 115:214 # Middle segment
  M2 <- 315:414 
  
  iSetup<-commSetup(S=S,L=L,W=W,compType="temp",years=iYears+ccYears,tempYSD=tempYSD,tempLSD=tempLSD)
  P<-iSetup$P   # Biotic parameters
  X<-iSetup$X   # Abiotic parameters
  Xi<-X # Save default X
  
  
  # Make a matrix for the variations on strategies
  J1 <- 9
  J2 <- 8
  J <- J1+J2-1
  jM <- matrix(c(c(seq(0,1,length=J1),seq(1,8,length=J2)[2:J2]),
                 c(seq(0,1,length=J1),seq(1,8,length=J2)[2:J2]),
                 c(seq(0,1,length=J1),seq(1,8,length=J2)[2:J2]),
                 c(seq(0,1,length=J1),seq(1,8,length=J2)[2:J2]),
                 c(seq(0,8,length=J1),seq(8,64,length=J2)[2:J2]),
                 c(seq(0,8,length=J1),seq(8,64,length=J2)[2:J2])),J,8)
  
  # Pre-allocate the matrix where everything is saved
  cOut <- matrix(NA,2*(1+6*(J-1)),39)
  sOut <- matrix(NA,2*(1+6*(J-1))*S,65)
  commList<-rep( list(list()), 5)

  k <- 0 # Keep track of the current row when saving the simulation
  l <- 0 # same for species row
  for(i in 1:2){
    # Do a simulation for each landscape
    q <- c(Q[,i]) # Current landscape
    X<-QChange(Xi,q,P$S) # Change to it
    Xii<-X # Save default X

    ni <- array(4,c(P$S,X$L,X$W)) # 4 individuals of all species in all locations... most won't survive
    iSim <- commSimulate(ni,P,X,y=1,years=iYears,init=T) # Initialize
    nif <- iSim$n # Initialized community

    # No action simulation
    k<-k+1
    nSim <- commSimulate(nif,P,X,y=iYears+1,years=ccYears)
    cOut[k,] <- commOut(id,tempLSD,tempYSD,iSim,nSim,X,P,2*i+1,0,0,M1,M2)
    for(s in 1:S){
      l<-l+1
      sOut[l,] <- specOut(id,s,tempLSD,tempYSD,iSim,nSim,X,P,2*i+1,0,0,M1,M2)
    }

    for(j in 2:J){
      # Corridor simulation
      k<-k+1
      Qc <- corridors(Xii$Q,jM[j,1])
      X<-QChange(Xii,Qc,S)
      cSim <- commSimulate(nif,P,X,y=iYears+1,years=ccYears)
      cOut[k,] <- commOut(id,tempLSD,tempYSD,iSim,cSim,X,P,2*i+1,1,jM[j,1],M1,M2)
      for(s in 1:S){
        l<-l+1
        sOut[l,] <- specOut(id,s,tempLSD,tempYSD,iSim,cSim,X,P,2*i+1,1,jM[j,1],M1,M2)
      }

      # Stepping stone simulation
      k<-k+1
      Qs <- steppingStone(Xii$Q,jM[j,2])
      X<-QChange(Xii,Qs,S)
      sSim <- commSimulate(nif,P,X,y=iYears+1,years=ccYears)
      cOut[k,] <- commOut(id,tempLSD,tempYSD,iSim,sSim,X,P,2*i+1,2,jM[j,2],M1,M2)
      for(s in 1:S){
        l<-l+1
        sOut[l,] <- specOut(id,s,tempLSD,tempYSD,iSim,sSim,X,P,2*i+1,2,jM[j,2],M1,M2)
      }

      # Reinforce simulation
      k<-k+1
      Qr <- reinforce(Xii$Q,jM[j,3])
      X<-QChange(Xii,Qr,S)
      rSim <- commSimulate(nif,P,X,y=iYears+1,years=ccYears)
      cOut[k,] <- commOut(id,tempLSD,tempYSD,iSim,rSim,X,P,2*i+1,3,jM[j,3],M1,M2)
      for(s in 1:S){
        l<-l+1
        sOut[l,] <- specOut(id,s,tempLSD,tempYSD,iSim,rSim,X,P,2*i+1,3,jM[j,3],M1,M2)
      }


      # All restore simulation
      k<-k+1
      Qa <- allRestore(Xii$Q,jM[j,4])
      X<-QChange(Xii,Qa,S)
      aSim <- commSimulate(nif,P,X,y=iYears+1,years=ccYears)
      cOut[k,] <- commOut(id,tempLSD,tempYSD,iSim,aSim,X,P,2*i+1,5,jM[j,4],M1,M2)
      for(s in 1:S){
        l<-l+1
        sOut[l,] <- specOut(id,s,tempLSD,tempYSD,iSim,aSim,X,P,2*i+1,5,jM[j,4],M1,M2)
      }


      # AM 50 simulation
      k<-k+1
      X<-Xii
      AM<-amSetup(P,X,ccYears,eta=50,AM=T,movesLeft=jM[j,5])
      amSim50<-commSimulate(nif,P,X,y=iYears+1,years=ccYears,AM=AM)
      cOut[k,] <- commOut(id,tempLSD,tempYSD,iSim,amSim50,X,P,2*i+1,7,jM[j,5],M1,M2)
      for(s in 1:S){
        l<-l+1
        sOut[l,] <- specOut(id,s,tempLSD,tempYSD,iSim,amSim50,X,P,2*i+1,7,jM[j,5],M1,M2)
      }


      # AM 75 simulation
      k<-k+1
      AM<-amSetup(P,X,ccYears,eta=75,AM=T,movesLeft=jM[j,6])
      amSim75<-commSimulate(nif,P,X,y=iYears+1,years=ccYears,AM=AM)
      cOut[k,] <- commOut(id,tempLSD,tempYSD,iSim,amSim75,X,P,2*i+1,8,jM[j,6],M1,M2)
      for(s in 1:S){
        l<-l+1
        sOut[l,] <- specOut(id,s,tempLSD,tempYSD,iSim,amSim75,X,P,2*i+1,8,jM[j,6],M1,M2)
      }
    }
    xs<-list(tempY=X$tempY,tempYSD=X$tempYSD,temp1d=X$temp1d,temp2d=X$temp2d)
    ps<-list(S=P$S,z=P$z,gam=P$gam,sig=P$sig,lam=P$lam,ro=P$ro,zo=P$zo)
    commList[[i]]<-list(X=xs,P=ps,iSim=iSim,nSim=nSim)

  }
  outList<-list(cOut=cOut,sOut=sOut,commList=commList)
  return(outList)
}


## below is code to run two iterations with several variations on a computer cluster
## to run a single simulation, use the simulation() function

# args<-commandArgs(TRUE)
# eval(parse(text=args[[1]]))
# 
# iF<-5000+i
# 
# cl<-makeCluster(4)
# funcList<-c(ls(),'erf','trapz','isempty','ceil','ifft')
# clusterExport(cl,funcList)
# outList <- clusterMap(cl,simulation,c(i,iF))
# stopCluster(cl)
# 
# cOuts<-outList[[1]]$cOut
# sOuts<-outList[[1]]$sOut
# 
#   cOuts<-rbind(cOuts,outList[[2]]$cOut)
#   sOuts<-rbind(sOuts,outList[[2]]$sOut)
# 
# write.table(cOuts,file=paste('corCommSept/comm',sprintf('%04d',i),".csv",sep=''),sep=',',col.names=F,row.names=F)
# write.table(sOuts,file=paste('corSpecSept/spec',sprintf('%04d',i),".csv",sep=''),sep=',',col.names=F,row.names=F)
# save(outList,file=paste('corOutSept/ol', sprintf('%04d',i),'.RData',sep=''))
