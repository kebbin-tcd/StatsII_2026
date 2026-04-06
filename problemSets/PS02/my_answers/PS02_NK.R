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

# here is where you load any necessary packages
# ex: stringr
# lapply(c("stringr"),  pkgTest)

lapply(c(),  pkgTest)

# set wd for current folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#####################
# Problem 1
#####################

# load data
load(url("https://github.com/ASDS-TCD/StatsII_2026/blob/main/datasets/climateSupport.RData?raw=true"))
#ls()
#class(climateSupport)
df <- as.data.frame(climateSupport)
summary(df)
#factoring variables
df$countries <- as.factor(df$countries)
df$sanctions <- as.factor(df$sanctions)
#question 1
predictive_model <- glm( #this is the predictive model
  choice ~ countries + sanctions,
  data = df,
  family = binomial(link = "logit")
)
summary(predictive_model)
null_model <- glm(choice ~ 1, data = df, family = binomial) #null model
#now to get the likelihood model we run anova
anova(null_model, predictive_model, test = "LRT")
#The pvalue is very small, smaller than significance level of 0.05, then we reject the null hypothesis. 
#The Deviance is 215.15, which is the improvement in fit when adding our variable to the null model.


#question 2
#part a
#new data frames
data_5_percent_160of192 <- data.frame(countries = factor("160 of 192", levels = levels(df$countries)), 
                                 sanctions = factor("5%", levels = levels(df$sanctions)))
data_15_percent_160of192 <- data.frame(countries = factor("160 of 192", levels = levels(df$countries)), 
                                  sanctions = factor("15%", levels = levels(df$sanctions)))
#prediction
log_odds_5_percent <- predict(predictive_model, newdata = data_5_percent_160of192, type = "link")
log_odds_15_percent <- predict(predictive_model, newdata = data_15_percent_160of192, type = "link")
change_in_log_odds <- log_odds_15_percent - log_odds_5_percent
odds_ratio <- exp(change_in_log_odds)
#The odds ratio for data section 5% to 15% sanctions at 160 countries is 0.722.
#An odds ratio > 1 means the odds increase.


#part b
#new data frames
data_5_percent_20of192 <- data.frame(countries = factor("20 of 192", levels = levels(df$countries)), 
                                      sanctions = factor("5%", levels = levels(df$sanctions)))
data_15_percent_20of192 <- data.frame(countries = factor("20 of 192", levels = levels(df$countries)), 
                                       sanctions = factor("15%", levels = levels(df$sanctions)))
#prediction
log_odds_5_percent_20 <- predict(predictive_model, newdata = data_5_percent_20of192, type = "link")
log_odds_15_percent_20 <- predict(predictive_model, newdata = data_15_percent_20of192, type = "link")
change_in_log_odds_20 <- log_odds_15_percent_20 - log_odds_5_percent_20
odds_ratio_20 <- exp(change_in_log_odds_20)
#The odds ratio for data section 5% to 15% sanctions at 20 countries is 0.722.
#An odds ratio > 1 means the odds increase.

#part c
#new data frames
data_80 <- data.frame(countries = factor("80 of 192", levels = levels(df$countries)), 
                               sanctions = factor("None", levels = levels(df$sanctions)))
#prediction
predicted_probability <- predict(predictive_model, newdata = data_80, type = "response")
#The estimate probability for 80 countries with no sanctions is 0.516

#question 3
#fitting a model with interaction
interaction_model <- glm(
  choice ~ countries * sanctions, #* for internaction
  data = df,
  family = binomial(link = "logit")
)
#comparing the two models
interaction_test <- anova(predictive_model, interaction_model, test = "LRT")
print(interaction_test)
#Pvalue is 0.39 which is higher than 0.05 so we fail to reject null hypothesis

