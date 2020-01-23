## SOURCE FUNCTIONS FOR KONA CRAB DEMOGRAPHIC ANALYSIS
##M KAPUR APR 2017
## path edits from Siple in Jan 2020

library(ggplot2)
library(Rmisc)
library(here)

## FUNCTIONS TO ESTIMATE DEMOGRAPHIC PARAMETERS TO INPUT TO LESLIE MATRIX

source(here("Maia code","sizeMatrix.R")) ## generate the size matrix separately

## FUNCTION ESTIMATE.AGE - INTERNAL USE ONLY
## PURPOSE: Generate a rough estimate of size at age for internal use only

estimate.age = function(size.bins, Linf, K, tZero, longevDraw){
  age.est.vec = NULL
  age.est.vec[ncol(size.bins)] = longevDraw ## the last value is always the longevity
  for (i in 1:(ncol(size.bins)-1)) {
    tau[i] = mean(size.bins[, i])
    j1[i] = min(size.bins[,i])+1E-10
    j2[i] = max(size.bins[,i])-1E-10
    age.est.vec[i] = (-1 / K) * log(1 - (j1[i] / Linf)) + tZero ## inverse VonB - this ensures it is self-dependent
  }
  # age.est.vec[ncol(size.bins)] = longevDraw ## the last value is always the longevity
  return(age.est.vec)
}
## FUNCTION HOENIG - INTERNAL USE ONLY
## PURPOSE: calculates mh for use in natural mortality
## O'Neill MH is 0.277; calculated from Hoening 1983 (lnZ = a + b lnTmax) with max longevity = 16y & using the "all"
# slope & intercept.
## For fish this value is 0.26.
hoenig = function(hoenig.slope,hoenig.int,longevDraw){
  mh = exp(hoenig.int + hoenig.slope*log(longevDraw))
  return(mh)
}
## FUNCTION NATURAL MORTALITY & FUNCTION SURVIVORSHIP - O'NEILL (2010, Spanner crabs)
## PURPOSE: generates natural mortality & survivorship at age
nat.mort = function(longevDraw, zeta, mh, age.est.vec) {
  M = NULL
  S = NULL
  for (a in 1:length(age.est.vec)) {
    M[a] = ((age.est.vec[a] + 1) ^ -zeta) / (longevDraw ^ -zeta) * mh
    S[a] = exp(-M[a])*0.064 ## ! larval mortality from Sakai 1971
  }
  return(data.frame('MORTALITY' = M, 'SURVIVORSHIP' = S))
}
## FUNCTION: FX
## PURPOSE: GENERATE EXPECTED EGGS @ SIZE BIN
## Equation given by Onizuka fecundity-size relationship, mature sizes only
FX.func = function(size.bins){
  eggs = NULL
  for(f in 1:ncol(size.bins)){
    if(mean(size.bins[,f]) > 60){
      eggs[f] = (beta*mean(size.bins[,f]) - 110433.1) ## in 1000s of Individuals
    } else {
      eggs[f] = 0 }
  }
  return(eggs)
}
## FUNCTION: OX
## PURPOSE: derive the proportion of mature females at size bin. Have to go through age first.
## APPROACH I: PER CARVALHO 2014
OX.funcI = function(L50, beta, size.bins, K, tZero, age.est.vec) {
  prop.mat = NULL
  tau = NULL
  prop.mat[ncol(size.bins)] = 1
  XM50 = (-1 / K) * log(1 - (L50 / Linf)) + tZero ## inverse vonb to get expected age at L50
  # XM50 = 5
  for (i in 1:(ncol(size.bins)-1)) {
    prop.mat[i] = 1 / (1 + exp(-beta * (age.est.vec[i] - XM50)))
  }
  return(prop.mat)
}
## FUNCTION: makeHarvestMat
## ! Reformatted from original to accept strategies. Make sure that the legal vec matches up.
makeHarvestMat = function(size.bins, harvConst, tc){
  harvestVec = rep(NA, ncol(size.bins))
  for(i in 1:ncol(size.bins)){
    if(max(size.bins[,i] < tc)){
      harvestVec[i] = 1 ## 100% survivorship for those smaller than selected for
    } else {
      harvestVec[i] = 1 - harvConst ## will apply whatever fishing pressure this is to all selected classes
    }
  }
  # Define harvest matrix
  HarvestMat = matrix(0,ncol(size.bins),ncol(size.bins))
  # Populate diagonal elements with "HarvestVec"
  diag(HarvestMat) = harvestVec
  return(HarvestMat)
}
## FUNCTION: runBissect ! NOT UPDATED
# Purpose: Run the BISSECTION to estimate the stationary harvest (see lecture 7 of FISH 458 by Andre Punt for BM)
runBissect = function(Ulow,Uhigh,bissectIters,bissectConv,LeslieMat,tc)
{
  bissectN = 0
  eigenDomCurr = 0
  while(bissectN < bissectIters & (abs(eigenDomCurr-1) > bissectConv))
  {
    bissectN = bissectN + 1
    # Compute Hcur
    Ucur = (Ulow+Uhigh)/2
    # Compute the eigenvalue for Ucur
    harvConstTemp = Ucur
    HarvestMat = makeHarvestMat(size.bins, harvConst = harvConstTemp, tc = tc) ## ASSUMES LEGAL ONLY
    AHmat = LeslieMat%*%HarvestMat
    eigenDomHcurr = max(abs(eigen(AHmat)$values))
    if (eigenDomHcurr < 1) (Uhigh = Ucur) ## pop is constant at eigen = 1
    if (eigenDomHcurr > 1) (Ulow = Ucur)
    Us = Ucur
    #if(Us==Uhigh) (Us = c(NA))
  }
  return(Us)
}
## FUNCTIONS TO SAVE AND PLOT MODELLED OUTPUTS
# FUNCTION: writeOutMC
# Purpose: save list "OutLeslieMC" to txt files
writeOutMC <- function(OutLeslieMC, Name)
{
  Path <- "outputs/"
  line = as.character(date()) ## include a timestamp
  ## special internal function for init list
  fnlist = function(x,fil){
    z = deparse(substitute(x))
    cat(unlist(z), "\n", file=fil)
    nams=names(x)
    for (i in seq_along(x) ){
      cat(nams[i], x[[i]], "\n",
          file=fil, append=TRUE)}
    write(paste0('NSIMS ',nsims.master), fil, append=TRUE)
  }
  fnlist(inits, paste(Path, Name, "_inits.txt", sep=''))
  ## OTHER INITS AND OUTPUTS
  write.table(OutLeslieMC$lifehist, paste(Path, Name, "_lifehist.txt", sep=''), sep=",", row.names=F) ## write lifhis
  write.table(OutLeslieMC$params, paste(Path, Name, "_params.txt", sep=''), sep=",", row.names=F)
  write.table(OutLeslieMC$outUs, paste(Path, Name, "_outUS.txt", sep=''), sep=",", row.names=F)
  ## GENERATE RISK TABLE
  ## Purpose: calculate proportion of negative R values for each combination of harvest & tc scenarios, write to
  table
  risktable = data.frame(matrix(NA, ncol = length(h.vector), nrow = length(tc.vector)))
  rownames(risktable) = paste(tc.vector)
  names(risktable) = h.vector
  for(t in 1:length(tc.vector)){
    for(h in 1:length(h.vector)){ ## loop through TCs and harvest breaks
      sub = subset(OutLeslieMC$params, harvConst == h.vector[h] & tc == tc.vector[t])
      # risktable[t,h] = mean(sub$rVal)
      risktable[t,h] = nrow(sub[sub$rVal < 0,])/nrow(sub)
    }
  }
  write.table(risktable, paste('outputs/', Name, "_risktable.txt", sep=''), sep=",", row.names=F) ##
  # write to file
}
# # FUNCTION: genRiskTable
# # Purpose: calculate proportion of negative R values for each combination of harvest & tc scenarios, write to table
# genRiskTable = function(OutLeslieMC, harvest.breaks, tc.vector){
# risktable = data.frame(matrix(NA, ncol = 10, nrow = length(tc.vector)))
# rownames(risktable) = paste(tc.vector)
# names(risktable) = seq(0,0.9,1/harvest.breaks)
# for(t in 1:length(tc.vector)){
# for(h in 1:length(seq(0,0.9,1/harvest.breaks))){ ## loop through TCs and harvest breaks
# sub = subset(test$params, harvConst == seq(0,0.9,1/harvest.breaks)[h] & tc == tc.vector[t])
# # risktable[t,h] = mean(sub$rVal)
# risktable[t,h] = nrow(sub[sub$rVal < 0,])/nrow(sub)
# }
# }
# write.table(risktable, paste('G:/Kona Crab/outputs/', Name, "_risktable.txt", sep=''), sep=",", row.names=F) ##
# write to file
# return(risktable)
# }
## PLOT EIGENVALUE - UNFISHED
makeHistEigenval = function(OutLeslieMC){
  ggplot(OutLeslieMC$params[OutLeslieMC$params['harvConst'] == 0,], aes(x = eigenDomAH))+
    theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    xlab('Eigenvalue')+
    ylab('Rel.Frequency')+
    # ggtitle('Eigenvalue') +
    scale_x_continuous(limits = c(0,4), breaks = seq(0,4,0.5))+
    geom_histogram(alpha = 0.5, bins = 50)+
    geom_vline(xintercept = 1, col = 'red', linetype = 'dotdash')+
    geom_vline(xintercept = median(OutLeslieMC$param[OutLeslieMC$params['harvConst'] == 0,]$eigenDomAH), linetype =
                 'dotted')
}
## PLOT TDOUBLE - UNFISHED
makeHistTdouble <- function(OutLeslieMC){
  ggplot(OutLeslieMC$params[OutLeslieMC$params['harvConst'] == 0,], aes(x = tDouble)) +
    theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    xlab('Population Doubling Time (years)')+
    ylab('')+
    # ggtitle('Population Doubling Time (years)') +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.5))+
    geom_histogram(alpha = 0.5, bins = 50)+
    geom_vline(xintercept = median(OutLeslieMC$params[OutLeslieMC$params['harvConst'] == 0,]$tDouble), linetype =
                 'dotted')
}
## PLOT RZERO - UNFISHED
makeHistRzero <- function(OutLeslieMC){
  ggplot(OutLeslieMC$params[OutLeslieMC$params['harvConst'] == 0,], aes(x = Rzero)) +
    theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    xlab('Net reproductive rate (R0)')+
    ylab('Density')+
    # ggtitle('Net reproductive rate (R0)') +
    # scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1))+
    geom_histogram(alpha = 0.5, bins = 50)+
    geom_vline(xintercept = median(OutLeslieMC$params[OutLeslieMC$params['harvConst'] == 0,]$Rzero), linetype =
                 'dotted')
}
## PLOT GENTIME - UNFISHED
makeHistGenTime <- function(OutLeslieMC){
  ggplot(OutLeslieMC$params[OutLeslieMC$params['harvConst'] == 0,], aes(x = genTime)) +
    theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    xlab('Generation time (years)')+
    ylab('')+
    # ggtitle('Generation time, T (years)') +
    scale_x_continuous(limits = c(3.01,3.05), breaks = seq(3.01,3.05,0.02))+
    geom_histogram(alpha = 0.5, bins = 50)+
    geom_vline(xintercept = median(OutLeslieMC$params[OutLeslieMC$params['harvConst'] == 0,]$genTime), linetype =
                 'dotted')
}
# ## PLOT RVAL
makeHistRval <- function(OutLeslieMC){
  ggplot(OutLeslieMC$params[OutLeslieMC$params['harvConst'] == 0,], aes(x = rVal)) +
    theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    xlab('Population Growth Rate (r)')+
    ylab('Density')+
    # ggtitle('Population Growth Rate (r), (Year ^-1)') +
    scale_x_continuous(limits = c(0,2), breaks = seq(0,2,0.5))+
    geom_histogram(alpha = 0.5, bins = 50)+
    # geom_histogram(data = comparison$params$rVal, alpha = 0.5, fill = 'blue', bins = 50) +
    geom_vline(xintercept = median(OutLeslieMC$params[OutLeslieMC$params['harvConst'] == 0,]$rVal), linetype =
                 'dotted')
}
## PLOT all of the above at once
param.plots = function(OutLeslieMC){
  multiplot(
    # makeHistEigenval(OutLeslieMC), ## dropped this cause R val is the same
    makeHistRval(OutLeslieMC),
    makeHistRzero(OutLeslieMC),
    makeHistTdouble(OutLeslieMC),
    makeHistGenTime(OutLeslieMC),
    cols = 2
  )
}
## DIAGNOSTIC PLOTS
# ELASTICITY ANALYSIS
right = function (string, char){
  substr(string,nchar(string)-(char-1),nchar(string))
}

overnam = function(OutLeslieMC){
  row.names(OutLeslieMC)[1:4] = right(row.names(OutLeslieMC[1:4,]),5) ## overwrite names
  row.names(OutLeslieMC)[5] = right(row.names(OutLeslieMC[5,]),6) ## overwrite names
  row.names(OutLeslieMC)[6] = right(row.names(OutLeslieMC[6,]),7) ## overwrite names
  row.names(OutLeslieMC)[7] = right(row.names(OutLeslieMC[7,]),4) ## overwrite names
  return(OutLeslieMC)
}
plotElast = function(OutLeslieMC, Name){
  StatsElastMC = as.data.frame(t(makeElastStats(OutLeslieMC, Name, write.file = F)))
  elastfec = StatsElastMC[1:7,]
  elastsurv = StatsElastMC[8:14,]
  elastfec = overnam(elastfec)
  elastsurv = overnam(elastsurv)
  fec =
    ggplot(NULL, aes()) +
    theme_bw()+
    ggtitle("Fecundity") +
    xlab('Size Bin (mm)') +
    ylab('') +
    scale_x_discrete(limits = c(row.names(elastfec)))+
    scale_y_continuous(limits = c(0,0.35), breaks = seq(0,0.35,0.05))+
    geom_bar(data = elastfec,stat = 'identity', colour = 'black', alpha = 0.5,
             aes(x = rownames(elastfec), y = elastfec[,'Mean'])) +
    theme( # remove the vertical grid lines
      panel.grid.major.x = element_blank() ,
      axis.text=element_text(size=12),
      axis.title=element_text(size=14),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line())
  surv =
    ggplot(NULL, aes()) +
    theme_bw()+
    ggtitle("Survivorship") +
    xlab('Size Bin (mm)') +
    ylab('') +
    scale_x_discrete(limits = c(row.names(elastsurv)))+
    scale_y_continuous(limits = c(0,0.35), breaks = seq(0,0.35,0.05))+
    geom_bar(data = elastsurv,stat = 'identity', colour = 'black', alpha = 0.5,
             aes(x = rownames(elastsurv), y = elastsurv[,'Mean'])) +
    theme( # remove the vertical grid lines
      panel.grid.major.x = element_blank() ,
      axis.text=element_text(size=12),
      axis.title=element_text(size=14),
      # explicitly set the horizontal lines (or they will disappear too)
      panel.grid.major.y = element_line())
  multiplot(fec,surv, cols = 2)
}
## plot the survivorship estimates over time (these are natural estimates and occur before H applied)
# plotSurv = function(OutLeslieMC, Name){
# as.data.frame(test$lifehist[,9:15]) %>% summarise_each(funs(mean))
#
# as.data.frame(t(test$lifehist[,9:15])) %>% group_by() %>% summarise_each(funs(mean))
#
# overnam(as.data.frame(t(test$lifehist[1:7,9:15])))
#
# means = NULL
# for(i in 9:15){
# means[i] = print(mean(test$lifehist[,i]))
# cbind(means, as.vector(row.names(overnam(as.data.frame(t(test$lifehist[1:7,9:15]))))))
# }
# ggplot(NULL, aes()) +
# theme_bw()+
# ggtitle("Survivorship") +
# xlab('Size Bin (mm)') +
# ylab('') +
# scale_x_discrete(limits = c(row.names(elastsurv)))+
# scale_y_continuous(limits = c(0,0.35), breaks = seq(0,0.35,0.05))+
# geom_bar(data = elastsurv,stat = 'identity', colour = 'black', alpha = 0.5,
# aes(x = rownames(elastsurv), y = elastsurv[,'Mean'])) +
# theme( # remove the vertical grid lines
# panel.grid.major.x = element_blank() ,
# axis.text=element_text(size=12),
# axis.title=element_text(size=14),
# # explicitly set the horizontal lines (or they will disappear too)
# panel.grid.major.y = element_line())
#
#
# }
## Save a plot of your choice to a file format of your choice; will be 8x10 inches
savePlots = function(plotfunc, OutLeslieMC, Name, vals, form, height, width){
  Path = "figures/"
  if(form == 'png'){
    png(paste0(Path, Name,"_",vals,'.',form), height = height, width = width)
  } else if(form == 'eps'){
    eps(paste0(Path, Name,"_",vals,'.',form), height = 8, width = 10)
  } else if(form == 'jpg'){
    jpeg(paste0(Path, Name,"_",vals,'.',form), height = 960, width = 1200)
  } else {
    print("form must be one of 'png','eps','jpg'")
  }
  plotfunc(OutLeslieMC)
  dev.off()
}
## function HARVESTPLOT
## purpose: plot misty mountains based on different H and TC
harvestPlot = function(OutLeslieMC.FILE, Name, form){
  Path = "figures/"
  if(form == 'png'){
    png(paste0(Path, Name,"_RvsHTC.",form), height = 960, width = 1200)
  } else if(form == 'eps'){
    eps(paste0(Path, Name,"_RvsHTC.",form), height = 8, width = 10)
  } else if(form == 'jpg'){
    jpeg(paste0(Path, Name,"_RvsHTC.",form), height = 960, width = 1200)
  } else {
    print("form must be one of 'png','eps','jpg'")
  }
  cols = c(grey(0.7, 0.5), grey(0.5, 0.5), grey(0.3, 0.5), grey(0.1, 0.5))
  if(is.character(OutLeslieMC.FILE)){
    params.df = read.table(OutLeslieMC.FILE, header = T,sep = ',')
  } else {
    params.df = OutLeslieMC.FILE
  }
  params.df$harvConst = as.factor(params.df$harvConst)
  par(mfrow = c(5, 2))
  par(oma = c(4, 4, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
  par(mar = c(1, 1, 2, 2)) # make the plots be closer together
  for (h in 1:length(h.vector)) {
    plot(1, type = "n",xlim = c(-1.5, 1.5),ylim = c(0, 8.5),frame.plot = F,
         main = paste0('Harvest Mortality = ', h.vector[h]), cex.main = 2,axes = F ) ## empty plot
    axis(side = 1,
         labels = (h == 9 | h == 10),
         cex.axis = 2) ## only bottom plots get x labels
    mtext('Population Growth Rate (r)',side = 1,outer = TRUE,line = 2,cex = 1.5) ## overall x title
    axis(side = 2,labels = (h %% 2 != 0),cex.axis = 2) ## only odd plots get y label
    mtext('Density',side = 2,outer = TRUE,line = 2,cex = 1.5)
    if (h == 1) {
      df0 = params.df[params.df$harvConst == h.vector[h], ]
      for (t in 1:length(tc.vector)) {
        df1 = df0[df0$tc == tc.vector[t], ]
        d = density(df1$rVal)
        polygon(d, col = cols[t], border = 'white')
        abline(v = 0, lwd = 2, lty = 2)
        lines(c(median(df1$rVal), median(df1$rVal)),
              c(0, 8),
              lwd = 1,
              lty = 2)
        # par(mai=c(0,0,0,0))
        legend("left",legend = paste0(unique(tc.vector), " mm"),text.width = max(sapply(text, strwidth)),
               col = cols,lwd = 5,cex = 2,horiz = F, bty = 'n', title = 'minimum capture size'
        )
      }
    } else{
      df0 = params.df[params.df$harvConst == h.vector[h], ]
      for (t in 1:length(tc.vector)) {
        df1 = df0[df0$tc == tc.vector[t], ]
        d = density(df1$rVal)
        polygon(d, col = cols[t], border = 'white')
        abline(v = 0, lwd = 2, lty = 2)
        lines(c(median(df1$rVal), median(df1$rVal)),c(0, 8),lwd = 1,lty = 2)
      }
    }
  }
}
## FUNCTION: makeParStats
## Purpose: Make table of descriptive statistics for demographic parameters & save to tile
makeParStats = function(OutLeslieMC, tc = 0, Name, write.file){
  OutLeslieMC.sub = OutLeslieMC$params[OutLeslieMC$params['harvConst'] == 0
                                       & OutLeslieMC$params['tc'] == tc ,] ## subset for UNFISHED
  # Define SUMMARY STATISTICS TABLE of MC simulations for demographic parameters
  NamesStats = c("Mean","Median","Lower 95% CI","Upper 95% CI","n")
  StatsParamsMC = data.frame(matrix(NA, ncol = ncol(OutLeslieMC.sub),
                                    nrow = length(NamesStats)))
  # Define variable names of summary statistics table
  names(StatsParamsMC) = names(OutLeslieMC.sub)
  row.names(StatsParamsMC) = NamesStats
  # Populate SUMMARY STATISTICS TABLE of MC simulations for demographic parameters
  for (i in 1:ncol(OutLeslieMC.sub)){
    meanVal = mean(OutLeslieMC.sub[,i])
    medianVal = median(OutLeslieMC.sub[,i])
    CI95L = quantile(OutLeslieMC.sub[,i],.025,names=FALSE)
    CI95U = quantile(OutLeslieMC.sub[,i],.975,names=FALSE)
    minVal = min(OutLeslieMC.sub[,i])
    maxVal = max(OutLeslieMC.sub[,i])
    n = signif(length(OutLeslieMC.sub[,i]),digits = 2)
    # Populate table
    StatsParamsMC[,i] = c(meanVal,medianVal,CI95L,CI95U,n)
  }
  if(write.file == T){
    Path <- "outputs/"
    write.table(StatsParamsMC, paste(Path, Name, "_StatsParamsMC.txt", sep=''), sep=",", row.names=F)
  }
  return(StatsParamsMC)
}
## FUNCTION: makeSADStats
## Purpose: Make table of descriptive statistics for SAD
makeSADStats = function(OutLeslieMC, tc, Name, write.file){
  OutLeslieMC.sub = OutLeslieMC$outSAD[OutLeslieMC$outSAD['harvConst'] == 0,] ## subset for UNFISHED
  # Define SUMMARY STATISTICS TABLE of MC simulations for SAD
  NamesStats = c("Mean","Median","Lower 95% CI","Upper 95% CI","n")
  StatsSADMC = data.frame(matrix(NA, ncol = ncol(OutLeslieMC.sub), nrow = length(NamesStats)))
  # Define variable names of summary statistics table
  names(StatsSADMC) = names(OutLeslieMC.sub)
  row.names(StatsSADMC) = NamesStats
  # Populate SUMMARY STATISTICS TABLE of MC simulations for SAD
  for (i in 1:ncol(OutLeslieMC.sub)){
    meanVal = mean(OutLeslieMC.sub[,i])
    medianVal = median(OutLeslieMC.sub[,i])
    CI95L = quantile(OutLeslieMC.sub[,i],.025,names=FALSE)
    CI95U = quantile(OutLeslieMC.sub[,i],.975,names=FALSE)
    minVal = min(OutLeslieMC.sub[,i])
    maxVal = max(OutLeslieMC.sub[,i])
    n = signif(length(OutLeslieMC.sub[,i]),digits = 2)
    # Populate table
    StatsSADMC[,i] = c(meanVal,medianVal,CI95L,CI95U,n)
  }
  if(write.file == T){
    Path <- "outputs/"
    write.table(StatsSADMC, paste(Path, Name, "_StatsSADMC.txt", sep=''), sep=",", row.names=F)
  }
  return(StatsSADMC)
}
## FUNCTION: makeoutElast.Stats
## Purpose: Make table of descriptive statistics for outElasticities
makeElastStats = function(OutLeslieMC, Name, write.file){
  OutLeslieMC.sub = OutLeslieMC$outElast[OutLeslieMC$outElast['harvConst'] == 0,] ## subset for UNFISHED
  # Define SUMMARY STATISTICS TABLE of MC simulations for elasticities
  NamesStats = c("Mean","Median","Lower 95% CI","Upper 95% CI","n")
  StatsElastMC = data.frame(matrix(NA, ncol = ncol(OutLeslieMC.sub), nrow = length(NamesStats)))
  # Define variable names of summary statistics table
  names(StatsElastMC) = names(OutLeslieMC.sub)
  row.names(StatsElastMC) = NamesStats
  # Populate SUMMARY STATISTICS TABLE of MC simulations for elasticities
  for(i in 1:ncol(OutLeslieMC.sub)){
    meanVal = mean(OutLeslieMC.sub[,i])
    medianVal = median(OutLeslieMC.sub[,i])
    CI95L = quantile(OutLeslieMC.sub[,i],.025,names=FALSE)
    CI95U = quantile(OutLeslieMC.sub[,i],.975,names=FALSE)
    minVal = min(OutLeslieMC.sub[,i])
    maxVal = max(OutLeslieMC.sub[,i])
    n = signif(length(OutLeslieMC.sub[,i]),digits = 2)
    # Populate table
    StatsElastMC[,i] = c(meanVal,medianVal,CI95L,CI95U,n)
  }
  if(write.file == T){
    Path <- "outputs/"
    write.table(StatsElastMC, paste(Path, Name, "_StatsElastMC.txt", sep=''), sep=",", row.names=F)
  }
  return(StatsElastMC)
}
makeElastStats = function(OutLeslieMC, Name, write.file){
  OutLeslieMC.sub = OutLeslieMC$outElast[OutLeslieMC$outElast['harvConst'] == 0,] ## subset for UNFISHED
  # Define SUMMARY STATISTICS TABLE of MC simulations for elasticities
  NamesStats = c("Mean","Median","Lower 95% CI","Upper 95% CI","n")
  StatsElastMC = data.frame(matrix(NA, ncol = ncol(OutLeslieMC.sub), nrow = length(NamesStats)))
  # Define variable names of summary statistics table
  names(StatsElastMC) = names(OutLeslieMC.sub)
  row.names(StatsElastMC) = NamesStats
  # Populate SUMMARY STATISTICS TABLE of MC simulations for elasticities
  for(i in 1:ncol(OutLeslieMC.sub)){
    meanVal = mean(OutLeslieMC.sub[,i])
    medianVal = median(OutLeslieMC.sub[,i])
    CI95L = quantile(OutLeslieMC.sub[,i],.025,names=FALSE)
    CI95U = quantile(OutLeslieMC.sub[,i],.975,names=FALSE)
    minVal = min(OutLeslieMC.sub[,i])
    maxVal = max(OutLeslieMC.sub[,i])
    n = signif(length(OutLeslieMC.sub[,i]),digits = 2)
    # Populate table
    StatsElastMC[,i] = c(meanVal,medianVal,CI95L,CI95U,n)
  }
  if(write.file == T){
    Path <- "outputs/"
    write.table(StatsElastMC, paste(Path, Name, "_StatsElastMC.txt", sep=''), sep=",", row.names=F)
  }
  return(StatsElastMC)
}
makeUsStats = function(OutLeslieMC, Name, write.file){
  tc.vector = tc.vector
  OutLeslieMC.sub = OutLeslieMC$outUs
  # Define SUMMARY STATISTICS TABLE of MC simulations for stable harvest rate
  NamesStats = c("Mean","Median","Lower 95% CI","Upper 95% CI",'min','max',"n")
  StatsUsMC = data.frame(matrix(NA, ncol = length(NamesStats), nrow = length(tc.vector)))
  # Define variable names of summary statistics table
  # names(StatsUsMC) = NamesStats
  row.names(StatsUsMC) = paste(tc.vector)
  # Populate SUMMARY STATISTICS TABLE
  StatsUsMC = data.frame(OutLeslieMC.sub %>% group_by(tc) %>% summarise('meanVal' = mean(Us),
                                                                        'medianVal' = median(Us),
                                                                        'CI95L' = quantile(Us,.025,names=FALSE),
                                                                        'CI95U' = quantile(Us,.975,names=FALSE),
                                                                        'minVal' = min(Us) ,
                                                                        'maxVal' = max(Us) ,
                                                                        'n' = n()))
  if(write.file == T){
    Path <- "outputs/"
    write.table(StatsUsMC, paste(Path, Name, "_StatsUsMC.txt", sep=''), sep=",", row.names=F)
  }
  return(StatsUsMC)
}