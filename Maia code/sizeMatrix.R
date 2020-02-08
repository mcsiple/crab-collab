## Generate size transition probability matrix following Siddeek et al (2016)
## Slight modification in calcluating Pmolt using empirical parameters. ! Indicates deviations.
# remove scientific numbers
options(scipen = 1)
#library(dplyr)

## This size.matrix will get multiplied by instantaneous survivorship to generate overall survivorship at size class
#source(here('Maia code','paramsScenarios.R'))

## This size.matrix will get multiplied by instantaneous survivorship to generate overall survivorship at size class.
#source('paramsScenarios.R')

## Generate size transition probability matrix following Siddeek et al (2016)
## 11mm is size @ first emergence, Brown (2008). The 58 -65 bin is the range from 0 to 80% fecundity, Onizuka (1972)
## Also Onizuka did not observe catches smaller than 35mm, but this covers the range
# size.bins = data.frame(matrix(c(11, 58, 59, 65, 66, 80, 81, 90, 91, 102, 103, Linf, Linf + 1, 150), nrow = 2))
# colnames(size.bins) = c('11-58', '59-65', '66-80', '81-90', '91-102','103-115','116+')

size.bins <- data.frame(matrix(c(11, 60, 61, 80, 81, 100, 101, 120, 121, 140,141, 160, 161, 180, 181, Linf, Linf + 1, 250), nrow = 2))
colnames(size.bins) = c('11-60','61-80', '81-100', '101-120', '121-140', '141-160','161-180','181-211','212+')

inv.logit <- function(x){
  #' @description takes the inverse logit of x
  #' @param x number to take inv logit of (num) 
  #' @return inverse logit of x, returns NA if is.na(x)
  rev <- exp(x)/(1+exp(x))
  return(rev)
}

## STEP 1: model individual variation in growth increment [Pij]
## empty matrix for probability densities of growth increment - "individual variability"
PIJ = matrix(NA, nrow = ncol(size.bins), ncol = ncol(size.bins))
## name bins
colnames(PIJ) = colnames(size.bins)
for (i in 1:nrow(PIJ)) {
  ## trap terminal value (115+)
  if (i == nrow(PIJ)) {
    PIJ[i, i] = 1 ## ! all terminal classes retained
  } else { ## prob of staying in size class
    ## tau is the midpoint of the contributing age class, which is taken from that col in size.bins
    tau = mean(size.bins[, i])
    ## mu is the expected growth increment based on that midpoint
    mu = achen + bchen * tau
    ## j1 is the min of the receiving age class
    j1 = min(size.bins[, i + 1])
    ## j2 is the max of the receiving age class
    j2 = max(size.bins[, i + 1])
    ## proportion that undergo growth increment
    PIJ[i, i] = 1 - pnorm(j2 - tau, mean = mu, sd = sdmolt)
    PIJ[i, i + 1] = pnorm(j2 - tau, mean = mu, sd = sdmolt) - pnorm(j1 - tau, mean = mu, sd = sdmolt)
  }
}#
# assign all lower triangle & upper values to 0 (it is impossible to shrink or skip a class)
PIJ[lower.tri(PIJ, diag = FALSE)] = 0
PIJ[is.na(PIJ)] = 0
## STEP 2: Generate vector of molt probatilities
m = NULL
for (i in 1:ncol(size.bins)) {
  tau[i] = mean(size.bins[, i])
  ## APPROACH I, as in Siddeek paper
  # m[i] = 1 / (1 + exp(csiddeek*(tau[i]-dsiddeek)))
  ## ! APPROACH II, find one-year molt probabilities for each size class based on Chen's estimate of 1-year molt
  #prob.
  m[i] = (0.0053*tau[i]^2 + 0.028*tau[i] + 6.35)^(-1) # MCS: changed this to intermolt period as function of size, based on a paper by Moksnes et al. 2014
 
  ## APPROACH III, all probabilities are I
  # m[i] = 1
}
## STEP 3: Generate transition matrix size.matrix XIJ
## generate empty matrix with column/row for each bin
size.matrix = matrix(NA, nrow = ncol(size.bins), ncol = ncol(size.bins))
colnames(size.matrix) = colnames(size.bins)
for (i in 1:ncol(size.bins)) {
  if (i == ncol(size.bins)) {
    ## trap terminal value
    size.matrix[i, i] = m[i] * PIJ[i, i] + (1 - m[i])
  } else {
    size.matrix[i, i] = m[i] * PIJ[i, i] + (1 - m[i])
    size.matrix[i, i + 1] = m[i] * PIJ[i, i + 1]
  }
}
size.matrix[is.na(size.matrix)] = 0
size.vector = NULL
for(i in 1:ncol(size.matrix)){
  if(i == ncol(size.matrix)){
    size.vector[i] = size.matrix[i,i-1]
  } else{
    size.vector[i] = size.matrix[i,i+1]
  }}
# ## VISUALIZATIONS
# par(mfrow = c(3, 1))
# par(mfrow = c(1, 1))
# ## FOR i = j and i,j both less than ncol
# ## make an empty plot with the right title
# xseq = seq(0, max(size.bins))
# plot(
# 1,
# type = "n",
# main = 'stay-within-class individual growth increment',
# ylim = c(0, 0.04),
# xlim = c(min(xseq), max(xseq))
# )
#
# for (k in 1:4) {
# ## loop through the first four rows/cols
# i = j = k ## i and j are equal
# tau = mean(size.bins[, i])## tau is the midpoint of the contributing age class
# ## mu is the expected growth increment based on that midpoint
# mu = achen + bchen * tau
# print(mu) ## this should vary
# ## j2 is the max of the receiving age class
# j2 = max(size.bins[, i + 1]) ## this should also vary
# ## generate the distribution appropriate for these params
# y = dnorm(xseq, mean = mu, sd = sdmolt)
# ## plot the curve given by these params on the graph
# lines(xseq, y, type = 'l')
# ## fill in the polygon based on the tau and j2
# xseq2 = seq(0, j2 - tau, length = 100)
# y2 = dnorm(xseq2, mean = mu, sd = sdmolt)
# polygon(c(0, xseq2, j2 - tau), c(0, y2, 0), col = rgb(1, 0, 0, 0.5))
# print(data.frame(
# 'MID' = tau,
# 'INC' = mu,
# 'PIJ' = pnorm(j2 - tau, mean = mu, sd = sdmolt)
# ))
# }
# abline(v = 0, lwd = 2)
#
# ## FOR N > J > I
# xseq = seq(0, max(size.bins))
# plot(
# 1,
# type = "n",
# main = 'next-class individual growth increment, i < j < n',
# ylim = c(0, 0.04),
# xlim = c(min(xseq), max(xseq))
# )
# j = 1
# for (i in 1:4) {
# for (j in 1:4) {
# if (i != j) {
# tau = mean(size.bins[, i])## tau is the midpoint of the contributing size class
# ## mu is the expected growth increment based on that midpoint
# mu = achen + bchen * tau
# print(mu) ## this should vary
# ## j1 is the min of the receiving age class
# j1 = min(size.bins[, i + 1])
# ## j2 is the max of the receiving age class
# j2 = max(size.bins[, i + 1]) ## this should also vary
# ## generate the distribution appropriate for these params
# y = dnorm(xseq, mean = mu, sd = sdmolt)
# ## plot the curve given by these params on the graph
# lines(xseq, y, type = 'l')
# ## fill in the polygon based on the tau and j2
# xseq2 = seq(j1 - tau, j2 - tau, length = 100)
# y2 = dnorm(xseq2, mean = mu, sd = sdmolt)
# polygon(c(j1 - tau, xseq2, j2 - tau), c(0, y2, 0), col = rgb(0, 1, 0, 0.5))
# print(data.frame(
# 'MID' = tau,
# 'J1' = j1,
# 'J2' = j2,
# 'INC' = mu,
# 'PIJ' = pnorm(j2 - tau, mean = mu, sd = sdmolt)
# ))
# } else if (i == j) {
# next
# }
# }
# }
# abline(v = 0, lwd = 2)
#
#
# ## FOR TERMINAL size
# xseq = seq(0, max(size.bins))
# plot(
# 1,
# type = "n",
# main = 'terminal-class individual growth increment',
# ylim = c(0, 0.04),
# xlim = c(min(xseq), max(xseq))
# )
# j = i = 5
# ## tau is the midpoint of the contributing age class, which is taken from that col in size.bins
# tau = mean(size.bins[, i])
# ## j1 is the min of the receiving age class
# j1 = min(size.bins[, i])
# ## mu is the expected growth increment based on that midpoint
# mu = achen + bchen * tau
# y = dnorm(xseq, mean = mu, sd = sdmolt)
# ## plot the curve given by these params on the graph
# lines(xseq, y, type = 'l')
# ## fill in the polygon based on the tau and j2
# xseq2 = seq(tau - j1, max(xseq), length = 100)
# y2 = dnorm(xseq2, mean = mu, sd = sdmolt)
# polygon(c(tau - j1, xseq2, max(xseq)), c(0, y2, 0), col = rgb(0, 0, 1, 0.5))
# print(data.frame(
# 'MID' = tau,
# 'INC' = mu,
# 'PIJ' = pnorm(tau - j1, mean = mu, sd = sdmolt)
# ))
# abline(v = 0, lwd = 2)