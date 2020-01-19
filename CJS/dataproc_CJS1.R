# PRocess data for CJS model
library(here)
library(tidyverse)
library(reshape2)
setwd(here("CJS"))

data <- read.csv(here("Heeia_crabdata_Apr2019.csv"))
str(data)

tags <- data %>% 
  filter(!is.na(Tag_Num)) %>%
  select(DeployDate, DeployTime, RetrieveDate,RetrieveTime,Sex,Keep_Release_Dead,Tag_YNNew,Tag_Num)

(nind <- length(unique(tags$Tag_Num)) )   # number of tagged individuals at start

yesses <- tags %>% 
  filter(Tag_YNNew == "Yes")

newtags <- tags %>%
  filter(Tag_YNNew == "New") 

# Individuals that were picked up but don't have tag records  
any(!yesses$Tag_Num %in% newtags$Tag_Num) 
addind <- which(!yesses$Tag_Num %in% newtags$Tag_Num) 
yesses$Tag_Num[addind]
newtags <- data.frame(DeployDate = NA,
                         DeployTime = NA,
                         RetrieveDate = "11-Mar-18", 
                         RetrieveTime = NA,
                         Sex = NA,
                         Keep_Release_Dead = "Release",
                         Tag_YNNew = "New",
                         Tag_Num = yesses$Tag_Num[addind])
newtags$Sex <- yesses$Sex[addind]

# Add fake entries to dataframe for tagged individuals
newtagdata <- bind_rows(tags,newtags)  %>% 
  mutate(time_per = NA)

# turn trap pickup dates into time periods for model - messy but works
datelookup <- data.frame(orig_date = unique(newtagdata$RetrieveDate),
                         time_per = c(1,1,1,2,3,3,4,4,5,5,6,6,7,7,8,9,1))

for (i in 1:nrow(newtagdata)){
  if(newtagdata$RetrieveDate[i] %in% datelookup$orig_date){
    pt <- which(datelookup$orig_date==newtagdata$RetrieveDate[i])
    newtagdata$time_per[i] <- datelookup$time_per[pt]
  }
}

unique(newtagdata$Tag_Num)
mt <- newtagdata %>%
  filter(Tag_YNNew =="Yes") %>%
  dcast(Tag_Num~time_per,fill=0) 

# First event: all tagged crabs that went in
newtagdata %>% filter(Tag_YNNew=="New")
dat <- newtagdata %>% 
  filter(Tag_YNNew=="New" & time_per == 1) %>% 
  select(Tag_Num,time_per) %>%
  mutate(t1 = 1) %>%
  merge(mt,all.x=TRUE)

 colnames(dat)[4:8] <- c("t2","t5","t6","t8","t9")
 dat <- dat %>% select(Tag_Num,t1,t2,t5,t6,t8,t9)
dat[is.na(dat)] <- 0
dat[-1,dat>1] <- 1
write.csv(dat,"CrabResight.csv")
