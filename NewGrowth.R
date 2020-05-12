# New growth parameter comparison
# Task: use data that were collected from a closer geographical location to HI. Originally we used parameters estimated in E Africa (Moksnes et al. 2015); Reviewer 3 recommended using Bonine et al. 2008 from Kosrae, Micronesia
library(here)
library(tidyverse)

dat <- read.csv(here('Heeia_crabdata_Apr2019.csv'))
head(dat)

# Check for length frequency bins
dat %>%
  filter(!Sex %in% c("FJ","J","NM")
         & !is.na(Sex)) %>%
  ggplot(aes(x=CWid_LR_mm)) + 
  geom_histogram(bins=30) +
  facet_wrap(~Sex) +
  theme_classic(base_size = 14) +
  xlab("Carapace width (mm)") +
  ylab("Count") +
  ggtitle("All crabs (2018-2019)")

max(dat$CWid_LR_mm,na.rm=T)

# HFP data do not have clear size frequency bins (maybe because of trap selectivity?) but we can use bins from Bonine et al. 2008. 


#  From Bonine et al 2008, max observed interval between release and recap (I in days) increased with carapace width (L in mm) according to I = exp(0.009+0.0314L)
Len <- 0:205
im <- data.frame(Len = Len,Min.Intermolt.Interval = exp(0.009+0.0314*Len))
plot(im$Len,im$Min.Intermolt.Interval, 
     ylab='Minimum intermolt interval (days)',
     xlab='Carapace width (mm)',
     xlim=c(0,200),pch = 19)

Len <- c( 104.3,122.6,143.5,162.4,177.6,189.0,203.1) # mean lengths of the size groupings for the age groups in Bonine et al. 2008 (SDs are also available in the paper)
# Time to first size class: 6 mos (182 days)
Time.At.Instar <- exp(0.009+0.0314*Len)


Age.Days <- vector()
#Age.Days[1] <- 0
Age.Days[1] <- 182
for(i in 2:length(Len)){
  Age.Days[i] <- Age.Days[i-1] + Time.At.Instar[i-1]
}


plot(Age.Days,Len,
     col='red',pch=19,
     type = 'p',
     xlim=c(0,1000),ylim=c(0,300),
     xlab = "Age (days)",
     ylab = "Carapace width (mm)")

t0 = -0.0015 # fix t0 
VBGM <- Len ~ Linf * (1 - exp(-K*(Age.Days-t0)))
VBGM_fit <- nls(VBGM, start = list (Linf = max(Len), K = 0.0057))
summary(VBGM_fit)

bonine.fits <- predict(VBGM_fit,newdata = list(Age.Days = 1:1500))
lines(1:1500,bonine.fits)

moksnes.ages <- 0:1500
moksnes.lengths <- 310*(1-exp(-0.57*((moksnes.ages/30)-(-0.019))))

lines(moksnes.ages,moksnes.lengths,col = 'blue')
legend('bottomright',pch = c(19,NA,NA),lty=c(0,1,1),legend = c('Size classes in Bonine et al. (2008)','Bonine et al. (2008) prediction','Moksnes et al. (2015) prediction'),col = c('red','black','blue'))


# Test with simulated data ------------------------------------------------
ages <- rep(seq(1,1100,50),each=10) #10 data pts at each age
detlengths <- 203.1*(1-exp(-0.0057*(ages-(-0.017))))
lengths <- rnorm(n=length(ages),mean = detlengths,sd=5)
plot(ages,lengths)
VBGM <- lengths ~ Linf * (1 - exp(-K*(ages-t0)))
VBGM_fit <- nls(VBGM, start = list (Linf = 203.1, K = 0.0057, t0 = -0.017))

summary(VBGM_fit)
# add_row(Len = 4, Time.At.Instar=1.144,Age.Days=1.144) # pre-molt is ~4 mm and min intermolt interval is 1.144 days



# Need to add some randomness ---------------------------------------------


