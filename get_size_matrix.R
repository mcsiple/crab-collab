# setup size matrix
dat <- read.csv(here::here('Heeia_crabdata_Apr2019.csv'))
(minsize <- min(dat$CWid_LR_mm, na.rm=T))
(maxsize <- max(dat$CWid_LR_mm, na.rm=T))

# Set up size bins
up_pts <- seq(minsize+4,maxsize,by=20)
mid_pts <- up_pts-10

paste0(mid_pts,up_pts)
# Generate size transition probability matrix following Siddeek et al (2016)