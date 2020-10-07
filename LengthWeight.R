# How to estimate length-weight parameters
library(here)
library(tidyverse)

allcrabs <- read.csv(here::here('Heeia_crabdata_Apr2019.csv'))

# Filter data to remove crabs that aren't mature enough to determine sex:
dat <- allcrabs %>%
  filter(!Sex %in% c("FJ","J","NM")
       & !is.na(Sex))

# Filter data to remove crabs that don't have a length or weight available:
dat <- dat %>%
  filter(!is.na(CWid_LR_mm) & !is.na(Weight_calc_kg))

# Calculate weight in g instead of kg
dat <- dat %>%
  mutate(Weight_calc_kg = as.numeric(Weight_calc_kg)) %>%
  mutate(Weight_calc_g = Weight_calc_kg * 1000)

# The classic equation for the relationship between length and weight is W = aL^b, where W is the weight and L is the "length" (in this case, cw).
# Technically this is W = aL^b * exp(epsilon), where epsilon is a normally distributed error term that describes the variation among individuals (you can think of W = aL^b  as being the average, with some individual variation around it). Because epsilon is normally distributed, exp(epsilon) is lognormally distributed.

# If you have an exponential relationship between x and y, you can log x and y and it is a linear relationship:
# These are just some made up values
a <- 0.0001
b <- 2.5

x <- 1:100
y <- a*x^b
plot(x,y)
plot(log(x),log(y))

# With some error:
x <- rep(1:100, each = 10)
for (i in 1:length(x)){
  epsilon <- rnorm(1,0,0.5)
  y[i] <- a * x[i] ^b * exp(epsilon)
}
plot(x,y)

plot(log(x),log(y))

# The easiest way to estimate the parameters a and b is to fit a linear model with lm(). 
# So plot log(length) on the x vs. log(weight) on the y for the crabs (I will leave this part for you to do!):



# To get the parameters that give you the best fit, you can use lm() with the log-transformed values:
log.x <- log(x)
log.y <- log(y)
fit <- lm(formula = log.y~log.x )
summary(fit) # in these estimates, the intercept is a and the slope is b
intercept = coef(fit)[1]
slope = coef(fit)[2]
predicted.log.y <- slope * log.x + intercept
lines(log.x,predicted.log.y,col='red')

# The estimates for a and b are still in log space so we have to transform them to get them back to normal.
plot(x,y)
a <- exp(intercept)
b <- slope
pred.y.normal <- a*x^b
lines(x,pred.y.normal,col='red')

# Voila! These are the a and b estimated from the data. Since we simulated the data in the first place, check to make sure a and b are the same as we hard-coded them on lines 25-26. 

# Now do the same thing for crabs, but with length as x and weight as y! Do it separately for males and females, because other studies have suggested their length-weight relationships are different. You can double check your estimates by plotting the fit on top of the data, like in lines 62-63.