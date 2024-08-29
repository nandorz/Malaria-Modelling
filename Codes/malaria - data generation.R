rm(list = ls())
# the whole code may take approximately 1.5 hours to run

################################################################################
# CHAPTER 1: DATA GENERATION FOR DAILY INCIDENCE ###############################
################################################################################

# 1. load data of yearly incidence from world malaria report Annex ####
library(pacman)
p_load(deSolve, tidyverse, gridExtra, readxl, ggplot2)

malaria.incidence <- read_excel("incidence_yearly.xlsx")
malaria.incidence <- as.data.frame(malaria.incidence)

# 2. taking the last 5 years data from the report ####
#pfpv.incidence <- malaria.incidence[1, c('2017', '2018', '2019', '2020', '2021')]
pf.incidence <- malaria.incidence[2, c('2017', '2018', '2019', '2020', '2021')]*0.95
#pv.incidence <- malaria.incidence[3, c('2017', '2018', '2019', '2020', '2021')]

# 3. preparing for turning yearly to daily incidence ####
# creating a multiplication factor for each day
# assuming the daily incidence has seasonality resembling cosine function
# for simplicity, the phase angle is set to be 0 so that the peak occurs at the beginning of the year
set.seed(123)
a <- c()
for(i in 1:360){
  a <- c(a,1+0.5*cos(2*pi*i/360) + rnorm(1, 0.1, 0.1))
}

# 4. make sure that the vector sum adds up to 1 ####
b <- a/sum(a)

# 5. multiplying the yearly incidence by daily multiplication factor ####
# pfpv.cases <- c(b*pfpv.incidence$'2017',b*pfpv.incidence$'2018', b*pfpv.incidence$'2019', b*pfpv.incidence$'2020', b*pfpv.incidence$'2021')
# setNames(pfpv.cases, 1:1800)

pf.cases <- c(b*pf.incidence$'2017',b*pf.incidence$'2018', b*pf.incidence$'2019', b*pf.incidence$'2020', b*pf.incidence$'2021')
setNames(pf.cases, 1:1800)

# pv.cases <- c(b*pv.incidence$'2017',b*pv.incidence$'2018', b*pv.incidence$'2019', b*pv.incidence$'2020', b*pv.incidence$'2021')
# setNames(pv.cases, 1:1800)

# 6. visual assessment of the daily incidence of the past 5 years ####
par(mfrow=c(1,1))
# plot(pfpv.cases, type="l", ylab = "population", xlab = "days", main = "P. falciparum and P. vivax")
plot(pf.cases, type="l", ylab = "population", xlab = "days", main = "P. falciparum")
# plot(pv.cases, type="l", ylab = "population", xlab = "days", main = "P. vivax")

# note that this data generation will be used for model fitting of susceptible compartment
# taking the whole population as the denominator would not make any sense as majority of indonesian area are not malaria-free
# however, the cases are clustered and concentration in certain areas, about 80% of cases in Papua region
# moreover, Indonesia consists of thousands of island so that cross-island transmission is less likely to be significant although not impossible
# For this reason, this model does not account for spatial analysis and imported cases
# Based on this generated data, a model fitting to estimate the number of susceptible population will be made
