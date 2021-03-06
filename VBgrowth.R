#Estimate some VB growth params

library(here)
library(reshape2)
library(tidyverse)
library(RColorBrewer)

#All sites together
all <- read.csv(here("Data","Heeia_crabdata_Apr2019.csv")



# Remove the one weird outlier
log.w <- log(all$Weight)
log.l <- log(all$Length)
lm1 <- lm(log.w~log.l)
plot(log.w~log.l,pch=19,xlab="Ln(length)",ylab="Ln(weight)",main="All Sites")
abline(lm1,col="red",lwd=2)

# Get a and b parameters from linear fit
int <- coef(lm1)[1]
slope <- coef(lm1)[2]
a <- exp(int)
b <- slope

# Try it using SSQ to see what happens (doesn't change results that much):
Lw.model <- function(params,length.vec=all$Length,weight.vec=all$Weight){
  a <- params[1]
  b <- params[2]
  pred.weights <- a * length.vec ^ b
  SSQ <- sum((weight.vec-pred.weights)^2)
  return(SSQ)
}

est <- optim(par=c(1,1),fn=Lw.model,method="Nelder-Mead")

# Fit length at age model 
all.growth <- subset(all, !is.na(Age)) # This takes out about 1000 entries that don't have age data
get.NLL <- function(params,Length.obs = all.growth$Length, age.vec = all.growth$Age) {
  Linf <- exp(params[1])
  k <- exp(params[2])
  t0 <- params[3]
  sigma <- params[4]
  L.pred <- Linf * (1-exp(-k*(age.vec-t0)))
  lt <- log(Length.obs)-log(L.pred)
  Lt <- (1/(sigma*sqrt(2*pi))) * exp( (-lt^2)/(2*sigma^2) )
  LL.vec <- log(Lt)
  NLL <- -1 * sum(LL.vec)
  return(NLL)
}

#Test starting values 
Linf = 300
k = 0.6
t0 = 1
sigma = 1

L.pred <- Linf * (1-exp(-k*(all.growth$Age-t0)))
#get.NLL(params = c(300,0.6,1,1))
optim(par=c(log(300),log(0.6),1,1),fn=get.NLL,method = "Nelder-Mead",control = list(maxit=100000))
# Est'd params
# Linf = 203.36
# k = 0.339
# t0 = -1.542
# sigma=0.078
par(mfrow=c(1,1))
plot(all.growth$Age,all.growth$Length,
     pch=19,col=rgb(0,0,0,alpha=0.3), 
     xlab="Age (years)",ylab="Length (mm)")
L.pred.plot <-203.36 * (1-exp(-0.339*(1:9+1.542)))  # Linf * (1-exp(-k*(all.growth$Age-t0)))
lines(1:9,L.pred.plot,col="red",lwd=2)

pred.table <- data.frame(age = 2:9, predlength = L.pred.plot[2:9] ) # Table of length predictions from age (for matching w growth data)
all.growth$preds <- pred.table[match(all.growth$Age,pred.table$age),'predlength']
all.growth$lengthresids <- all.growth$Length - all.growth$preds
ggplot(all.growth,aes(x=Age,y=lengthresids,colour=Site)) +
  geom_point(alpha=0.6,size=2) +
  scale_colour_brewer(type="qual",palette=3) +
  facet_wrap(~Site) +
  geom_hline(yintercept = 0,lty = 2) +
  ylab("Residuals from age-length relationship (mm)") +
  theme_bw(base_size = 14)

yearpalette <- brewer.pal(n=8,name = "RdYlGn")
ggplot(all.growth,aes(x=Age,y=lengthresids,colour=Year)) +
  geom_point(alpha=0.6,size=2) +
  scale_colour_gradientn(colours = yearpalette) +
  facet_wrap(~Year) +
  geom_hline(yintercept = 0,lty = 2) +
  ylab("Residuals from age-length relationship (mm)") +
  theme_bw(base_size = 14)

# Now get predicted weights at age from these fits, and add to prediction table:
pred.table$predweight <- a*pred.table$predlength^b
save(pred.table,file = "PredLengthsWeights_mm_g.RData")
# Check to see if weight at age is realistic, so you can use the weights at age to match

all.growth$weightpreds <- pred.table[match(all.growth$Age,pred.table$age),'predweight']
all.growth$weightresids <- all.growth$Weight-all.growth$weightpreds
ggplot(all.growth,aes(x=Age,y=weightresids,colour=Site)) +
  geom_point(alpha=0.6,size=2) +
  scale_colour_brewer(type="qual",palette=3) +
  facet_wrap(~Site) +
  geom_hline(yintercept = 0,lty = 2) +
  ylab("Residuals from weight at age relationship (g)") +
  theme_bw(base_size = 14)

