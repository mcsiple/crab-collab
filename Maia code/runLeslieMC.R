#runLeslieMC


runLeslieMC <- function(nsims.master, harvest.breaks) {
  # Define output table of MC simulations for demographic parameters
  params.to.save = c('eigenDomAH',
                     'rVal',
                     'tDouble',
                     'Rzero',
                     'genTime',
                     'harvConst',
                     'tc')
  OutParams = data.frame(matrix(
    NA,
    nrow = nsims.master * harvest.breaks * length(tc.vector),
    ncol = length(params.to.save)
  ))
  names(OutParams) = params.to.save
  
  NamesLifeHist <- names(size.bins)
  InLifeHistMC = matrix(NA,
                        nrow = nsims.master,
                        ncol = length(NamesLifeHist) * 2 + 1,
                        dimnames = list(NULL, c(
                          "simID",
                          paste('Fec', NamesLifeHist),
                          paste('Surv', NamesLifeHist)
                        ))
  )
  
  
  # Define output table of MC simulations for SAD
  OutSADMC = data.frame(matrix(
    NA,
    nrow = nsims.master * harvest.breaks,
    ncol = ncol(size.bins) + 1
  ))
  names(OutSADMC) = c(paste0(colnames(size.bins)), 'harvConst')
  
  # Define output table of MC simulations for elasticities
  NamesElast  = c(paste0('elastFec', colnames(size.bins)),
                  paste0('elastSurv', colnames(size.bins)),
                  'elastCheck')
  OutElastMC = data.frame(matrix(
    NA,
    nrow = nsims.master * length(h.vector) * length(tc.vector),
    ncol = length(NamesElast) + 1
  ))
  names(OutElastMC) =  c(paste0(NamesElast), 'harvConst')
  
  # Define output table of MC simulations for Us (sustainable harvest rate)
  OutUsMC <-
    data.frame(matrix(NA, nrow = nsims.master * length(tc.vector), 3))
  names(OutUsMC) = c('simID', 'Us', 'tc')
  
  # MONTE CARLO simulations
  for (isim in 1:nsims.master) {
    if (isim %% 100 == 0) {
      cat("isim ",isim,"\n")
    } ## report every 100 sims
    
    ## Run Pop Model
    longevDraw = runif(1,7,16) ## lit reported range for Kona crab - don't have values for Samoan crab
    # longevDraw = 16
    
    mh = hoenig(
      hoenig.slope = hoenig.slope,
      hoenig.int = hoenig.int,
      longevDraw = longevDraw
    )
    
    nat.survivorship.temp = nat.mort(
      longevDraw = longevDraw,
      zeta = 0.3,
      mh =mh,
      age.est.vec = 1:longevDraw
    )[, 'SURVIVORSHIP'] 
    
    
    # Natural Mortality: Randomly choose survivorship from uniform distribution, dims of longevity (which varies)
    #SurvAge <- sort(runif(longevDraw, 0.01, 0.04), decreasing = F)
    SurvAge <- nat.survivorship.temp
    #cat("SurvAge  ",SurvAge,"\n")
    # survPlot(SurvAge, size.matrix)
    
    ## N is placeholder array of 1
    N <- array(1, dim = c(length(SurvAge), ncol(size.matrix)))
    ## Multiply by STM - this becomes diagonal of Lefkovich
    
    surv.vector <- (N * SurvAge) %*% size.matrix %>%
      round(digits = 4) %>%
      apply(.,2,mean) %>%
      t() %>%
      as.data.frame()
    surv.vector[1] <- 0.0001 # from Quinitio et al 2001; daily survival .266--> annual survival 0.0001 (or even less - but have to change it )  (spanner crab value 0.064)
    colnames(surv.vector) <- NULL
    #print(surv.vector)
    ## Fecundity
    FX = FX.func(size.bins,size.at.maturity) ## just expected egg output
    # EPR = ceiling(abs(rnorm(1,1,1)))
    # cat(EPR,"\n")
    EPR = ceiling(abs(runif(1,0,3))) #2
    fecundity.vector = FX * SR * EPR ## per capita female birth rate
    # cat(fecundity.vector,"\n")
    
    
    # surv.vector = nat.survivorship.temp %*% size.matrix ## multiply by size matrix to get surviving portion within each size class; see Punt 2003. These are decoupled from harvest
    
    ## GENERATE VECTORS FOR ESTIMATED RECRUITS
    # %>% apply(., 1, prop.table) %>% round(digits = 4)
    #### LESLIE MATRIX ####
    LeslieMat = matrix(0,
                       nrow = length(surv.vector),
                       ncol = length(surv.vector))
    # Populate first row with fertility vector
    LeslieMat[1, ] = fecundity.vector
    # Populate the rest of the matrix with survival terms (Pi) - only on diagonal
    for (i in 2:length(surv.vector)) {
      j = i - 1
      LeslieMat[i, j] = surv.vector[,i - 1]
    }
    
    # Populate LH list with draws for each run
    InLifeHistMC[isim,] = c(isim,fecundity.vector,as.numeric(surv.vector))
    
    
    ## Estimate the STATIONARY HARVEST (Us) using the BISSECTION METHOD
    ## [only relies on Leslie Mat, which already has nat mort in place]
    for (t in 1:length(tc.vector)) {
      ## all ages, only L50, only legal size
      indexUs = (isim - 1) * length(tc.vector) + t  ## for OutUS
      
      Us = runBissect(Ulow,
                      Uhigh,
                      bissectIters,
                      bissectConv,
                      LeslieMat,
                      tc = tc.vector[t]) ## generate different Us estimates for each TC value for each sim of Leslie
      OutUsMC[indexUs,] = as.vector(c(isim, Us, tc.vector[t]))
    }
    
    
    ## LEGAL ONLY BISSECTION
    ## Estimate the STATIONARY HARVEST (Us) using the BISSECTION METHOD [only relies on Leslie Mat, indifferent to harvest rate]
    # Us = runBissect(Ulow,
    #                 Uhigh,
    #                 bissectIters,
    #                 bissectConv,
    #                 LeslieMat,
    #                 tc = tc.vector[t]) ## generate different Us estimates for each TC value for each sim of Leslie
    # OutUsMC[isim, ] = as.vector(c(isim, Us, tc.vector[t]))
    
    
    
    #### HARVEST MATRIX ####
    ## Use Decrement -- we want to apply the entire suite of harvest scenarios to a single leslie
    
    for (h in 1:length(h.vector)) {
      for (t in 1:length(tc.vector)) {
        ## all ages, only L50, only legal size
        # t = 4 ## hard assignment of tc vector for legal sizes
        HarvestMat =  makeHarvestMat(size.bins, harvConst = h.vector[h], tc = tc.vector[t])
        # index = (isim-1)*length(h.vector) + h ## original just with H
        index = (isim - 1) * length(h.vector) * length(tc.vector) + (h - 1) *
          length(tc.vector) + t ## for bigger MCs
        
        #### DEMOGRAPHIC PARAMETERS ####
        # Dominant EIGENVALUE of the AH matrix
        AHmat = LeslieMat %*% HarvestMat ## MATRIX MULTIPLICATION OPERATOR
        # cat(diag(AHmat),"\n")
        eigenDomAH = max(abs(eigen(AHmat)$values)) ## a value of zero means it is NOT invertible...
        
        # Instantaneous rate of population increase (r) [exp(r) = lambda]
        rVal = log(eigenDomAH)
        # if (h > 8) {cat(rVal,"\n")}
        # Population doubling time (t2) and halving time (t0.5)
        tDouble = log(2) / rVal
        tHalf = log(.5) / rVal
        
        # NET REPRODUCTIVE RATE, R0 (see Caswell 2001, p. 126)
        ## Decompose the AH matrix (see Caswell 2001, p. 110, eq. 5.1)
        # Compute the T matrix (describing transitions)
        Tmat = AHmat
        Tmat[1,] = 0
        # Compute the F matrix (describing reproduction)
        Fmat = matrix(0, ncol(LeslieMat), ncol(LeslieMat))
        Fmat[1,] = AHmat[1,]
        ## Compute the FUNDAMENTAL MATRIX (N) (see Caswell 2001, p. 112, eq. 5.7)
        # Compute the IDENTITY MATRIX (I)
        Imat = matrix(0, ncol(LeslieMat), ncol(LeslieMat))
        diag(Imat) = 1
        # Compute N matrix
        Nmat = solve(Imat - Tmat)
        ## Compute R MATRIX (R) (see Caswell 2001, p. 126, eq. 5.64)
        Rmat = Fmat %*% Nmat
        ## Compute NET REPRODUCTIVE RATE, R0 (see Caswell 2001, p. 126)
        Rzero = max(abs(eigen(Rmat)$values))
        
        # GENERATION TIME, Abar (see Caswell 2001, p. 128)
        # How to compute Abar? Use T meanwhile...
        genTime = log(Rzero) / log(eigenDomAH)
        
        # Populate vector with output demographic parameters
        # OutParams[isim, ] = c(eigenDomAH, rVal, tDouble, Rzero, genTime, h)
        OutParams[index, ] = c(eigenDomAH,
                               rVal,
                               tDouble,
                               Rzero,
                               genTime,
                               h.vector[h],
                               tc.vector[t])
        
        
        ### FIX PROBLEM - hitting the upper bound for high tc values!
        
        ## ELASTICITY ANALYSIS (see Caswell 2001, p. 128)
        # Compute STABLE AGE DISTRIBUTION, SAD (w)
        w = abs(eigen(AHmat)$vectors[, 1])
        # Normalize w
        w = w / sum(w)
        # Compute REPRODUCTIVE VALUE vector (v)
        v = abs(eigen(t(AHmat))$vectors[, 1])
        # Compute inner product <w,v>
        iProd = sum(w * v)
        
        ## Compute SENSITIVITY matrix (see Caswell 2001, p. 209, eq. 9.12)
        # Define dimensions of sensitivity matrix
        SensMat = matrix(0, ncol(size.bins), ncol(size.bins))
        # Populate matrix with sensitivity elements
        for (i in 1:ncol(size.bins)) {
          for (j in 1:ncol(size.bins))
            SensMat[i, j] = (v[i] * w[j]) / iProd
        }
        
        ## Compute ELASTICITY matrix (see Caswell 2001, p. 226, eq. 9.70)
        # Define dimensions of elasticity matrix
        ElastMat = matrix(0, ncol(size.bins), ncol(size.bins))
        # Populate matrix with sensitivity elements
        for (i in 1:ncol(size.bins)) {
          for (j in 1:ncol(size.bins))
            ##! this is wrt lambda - does it matter?
            ElastMat[i, j] = (AHmat[i, j] / eigenDomAH) * SensMat[i, j]
        }
        ## Compute FECUNDITY and SURVIVAL ELASTICITIES for age-groups
        # Extract first row (fecundity elasticities) from "ElastMat"
        ElastFec = ElastMat[1,]
        # Extract vector with survival elasticities from "ElastMat"
        ElastSurv = rep(0, ncol(size.bins))
        for (i in 2:ncol(size.bins)) {
          j = i - 1
          ElastSurv[i - 1] = ElastMat[i, j]
        }
        # meaningless names just placeholders
        elastFec_1 = ElastFec[1]
        elastFec_2 = ElastFec[2]
        elastFec_3 = ElastFec[3]
        elastFec_4 = ElastFec[4]
        elastFec_5 = ElastFec[5]
        elastFec_6 = ElastFec[6]
        elastFec_7 = ElastFec[7]
        elastFec_8 = ElastFec[8]
        elastFec_9 = ElastFec[9]
        
        elastSurv_1 = ElastSurv[1]
        elastSurv_2 = ElastSurv[2]
        elastSurv_3 = ElastSurv[3]
        elastSurv_4 = ElastSurv[4]
        elastSurv_5 = ElastSurv[5]
        elastSurv_6 = ElastSurv[6]
        elastSurv_7 = ElastSurv[7]
        elastSurv_8 = ElastSurv[8]
        elastSurv_9 = ElastSurv[9]
        
        elast.vals = c(
          elastFec_1,
          elastFec_2,
          elastFec_3,
          elastFec_4,
          elastFec_5,
          elastFec_6,
          elastFec_7,
          elastFec_8,
          elastFec_9,
          elastSurv_1,
          elastSurv_2,
          elastSurv_3,
          elastSurv_4,
          elastSurv_5,
          elastSurv_6,
          elastSurv_7,
          elastSurv_8,
          elastSurv_9
        )
        # Check on elasticity computations
        elastCheck = sum(ElastFec) + sum(ElastSurv)
        
        ## Populate output table of elasticity analysis
        OutElast = c(elast.vals, elastCheck)
        
        ## Populate output tables
        # Table of SAD
        OutSAD = rep(NA, length(OutSADMC))
        OutSAD = w
        OutSAD[ncol(size.bins)] = sum(w[ncol(size.bins) - 1:ncol(size.bins)])
        
        OutSADMC[index,] = as.vector(c(OutSAD, h.vector[h]))
        OutElastMC[index,] = as.vector(c(OutElast, h.vector[h]))
      }
    }
  }
  
  return(
    list(
      'params' = OutParams,
      'lifehist' = InLifeHistMC,
      'inits' = inits,
      'outSAD' = OutSADMC,
      'outElast' = OutElastMC,
      'outUs' = OutUsMC
    )
  )
}
