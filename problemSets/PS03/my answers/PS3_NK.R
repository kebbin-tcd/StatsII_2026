#####################
# load libraries
# set wd
# clear global .envir
#####################

# remove objects
rm(list=ls())
# detach all libraries
detachAllPackages <- function() {
  basic.packages <- c("package:stats", "package:graphics", "package:grDevices", "package:utils", "package:datasets", "package:methods", "package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:", search()))==1, TRUE, FALSE)]
  package.list <- setdiff(package.list, basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package,  character.only=TRUE)
}
detachAllPackages()

# load libraries
pkgTest <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[,  "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg,  dependencies = TRUE)
  sapply(pkg,  require,  character.only = TRUE)
}
library(dplyr)
# here is where you load any necessary packages
# ex: stringr
# lapply(c("stringr"),  pkgTest)

lapply(c("nnet", "MASS"),  pkgTest)
library(nnet)

# set wd for current folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#####################
# Problem 1
#####################

# load data
gdp_data <- read.csv("https://raw.githubusercontent.com/ASDS-TCD/StatsII_2026/main/datasets/gdpChange.csv", stringsAsFactors = F)
#part 1
head(gdp_data)
gdp_data$REG <- factor(gdp_data$REG,
                           levels = c("0","1"),
                           labels = c("Non-Democracy", 
                                      "Democracy"))
head(gdp_data)
gdp_data <- gdp_data %>%
  mutate(GDPWdiff_category = case_when(
    GDPWdiff < 0  ~ "negative",
    GDPWdiff == 0 ~ "no change",
    GDPWdiff > 0  ~ "positive"))
head(gdp_data)
gdp_data$GDPWdiff_category <- factor(gdp_data$GDPWdiff_category)
gdp_data$GDPWdiff_category <- relevel(gdp_data$GDPWdiff_category, ref = "no change")
multi.log <- multinom(GDPWdiff_category ~ REG + OIL, data = gdp_data)
summary(multi.log)
# run model
exp(coef(multi.log))
# get p values
z <- summary(multi.log)$coefficients
z
#interpretation
#Holding other factors constant, the odds of a country having positive GDP growth rather than no change are 5.86 times higher for democracy than for non-democrac.
#Holding other factors constant, the odds of a country having negative GDP growth rather than no change are 3.97 times higher for democracy
#The odds of positive GDP growth versus no change is 97.15 times higher for countries exporting oil
#The odds of negative GDP growth versus no change is 119.57 times higher for countries exporting oil

#part 2
gdp_data$GDPWdiff_category <- factor(gdp_data$GDPWdiff_category, 
                                     levels = c("negative", "no change", "positive"), 
                                     ordered = TRUE)
ord.log <- polr(GDPWdiff_category ~ REG + OIL, data = gdp_data, Hess = TRUE)
summary(ord.log)
exp(coef(ord.log))
#intpretation
#Holding other variables constant, the odds of a country being in a higher GDP growth are 1.48 times higher for democracrativ countries. The t-value of democracy is 5.3 hence highly significant.
#Holding other variables constant, the odds of a country being in a higher GDP growth are 0.8 times lowers for countries exporting oil. THe t-value of oil is -1.7 so it is not significant.

#####################
# Problem 2
#####################

# load data
mexico_elections <- read.csv("https://raw.githubusercontent.com/ASDS-TCD/StatsII_2026/main/datasets/MexicoMuniData.csv")
head(mexico_elections)
#part a
mexico_elections_poisson <- glm(PAN.visits.06 ~ competitive.district + marginality.06 + PAN.governor.06, 
                 data = mexico_elections, 
                 family = "poisson")
summary(mexico_elections_poisson)
#interpret
#The p-value for competitive district is 0.63, hence not very significant.
#part b
exp(coef(mexico_elections_poisson))
#For every 1 unit increase in the poverty measure, the expected number of visits decreases by by 88% holding other variables constant.
#For PAN affiliated governers they recieve 27% visits less than non-affiliated PAN governors holding other variables constant.

#pact c
hypothetical_district <- data.frame(competitive.district = 1, 
                                marginality.06 = 0, 
                                PAN.governor.06 = 1)
predict(mexico_elections_poisson, newdata = hypothetical_district, type = "response")
#The result is 0.01494818 