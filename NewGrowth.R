# New params attempt
library(here)
library(tidyverse)

dat <- read.csv(here('Heeia_crabdata_Apr2019.csv'))
head(dat)

# According to Christine and Kiva, you can get the VBGF parameters by plotting carapace width vs. summed intermolt time at each instar (i.e., )

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

#  From Bonine et al 2008, max observed interval between release and recap (I in days) increased with carapace width (L in mm) according to I = exp(0.009+0.0314L)
Len = c(104.3,122.6,143.5,162.4,177.6,189.0,203.1) # mean lengths of the size groupings (SDs are also available in the paper)
im <- data.frame(Len = Len,Interval = exp(0.009+0.0314*Len))
plot(im$Interval, im$Len,
     xlab='Maximum observed interval between release and recapture (days)',
     ylab='Carapace width (mm)')
