# analysis using repeated augmented dickey fuller test -----------------------------

# load useful libraries 
library(tseries)
library(xtable)
library(stats)
library(MASS)
library(dplyr, warn.conflicts= FALSE)
library(ggplot2)
library(GlmSimulatoR)
library(tidyverse)
library(ggfortify)
library(ggseg)
library(ggseg3d)


# load data:
load("data/fMRI-ROI-time-series.RData") 
covariates <- read.table("data/subj-covariates.txt", header = T)
covariates$type_current_diagnosis[is.na(covariates$type_current_diagnosis)] <- "-"
covariates$type_lifetime_diagnosis[is.na(covariates$type_lifetime_diagnosis)] <- "-"

# data cleaning
Y1 <- Y[,1:403,c(-1,-21),1]
covariates <- covariates[c(-1,-21),]

# run repeated DF test 

repeated_df <- function(Matrix){
  s<- NULL
  test <- NULL
  adf<- matrix(NA, nrow = 22, ncol = 70)
  for(i in 1: dim(Matrix)[3]){
    for(j in 1:dim(Matrix)[1]){
      s <- as.vector(na.omit(Matrix[j,,i]))
      test <- adf.test(s)
      adf[i,j] <- ifelse(test$p.value <= 0.01,1,0)
    }
  }
  n.stationary <- cbind(colnames(Matrix[,1,]), apply(adf, 1, sum))
  n.stationary # these numbers rapresent the number of (non) stationary series per subject
  rownames(n.stationary) <- n.stationary[,1]
  n.stationary <- as.data.frame(n.stationary[,-1])
  n.stationary.1 <- n.stationary[1:11,]
  n.stationary.2 <- n.stationary[12:22,]
  stationary.table <- cbind(rownames(n.stationary)[1:11], n.stationary.1, 
                            rownames(n.stationary)[12:22], n.stationary.2)
  rownames(stationary.table) <- NULL
  return(list(stationary.table, n.stationary))
}

results  <- repeated_df(Y1)
stationary.table <- results[[1]]
n.stationary <- results[[2]]

# export this table for latex
print(xtable(stationary.table, type = "latex"), file = "latex_elements/stationarytable.tex")



# Visually check the validity of results  by looking at the 24 different individuals 
# and see whether this statistic really shows a substantial and data analytically 
# clear difference between the time series 


# let's explore the situation for subject_9 (with n=62) against subject_10 (with n=70)

par(mfrow = c(1,3))
ts.plot(Y1[1,,8], ylim = c(-20,20), ylab = "BOLD")
for(i in 2:70){
  lines(Y1[i,,8], col=i)
} #subject 9
ts.plot(Y1[1,,9], ylim = c(-20,20), ylab = "BOLD")
for(i in 2:70){
  lines(Y1[i,,9], col=i)
}#subject 10
ts.plot(Y1[1,,12], ylim = c(-20,20), ylab = "BOLD")
for(i in 2:70){
  lines(Y1[i,,12], col=i)
}


# For a single individual I can even look at 5 series that look stationary and 5
# that don't look stationary, and ask myself do that really properly look 
# different and does an individual that has 15 of them look properly different 
# from an individual that has 2 of them. And if this is the case, maybe it makes 
# some sense. 

# Let's take subject 13

# brain subjs 4, 42, 45, 51, 52, 54 are stationary according to the adf test

par(mfrow=c(2, 4))
ts.plot(Y1[4,,12], ylab = "STATIONARY",
        main = "REGION 4",
        ylim = c(-20,25))
ts.plot(Y1[42,,12], main="REGION 42",
        ylim = c(-20,25))
ts.plot(Y1[45,,12], main = "REGION 45",
        ylim = c(-20,25))
ts.plot(Y1[51,,12], main = "REGION 51",
        ylim = c(-20,25))
ts.plot(Y1[3,,12],ylab = "NON-STATIONARY", main = "REGION 3",
        ylim = c(-20,25))
ts.plot(Y1[28,,12], main = "REGION 28",
        ylim = c(-20,25))
ts.plot(Y1[46,,12], main= "REGION 46",
        ylim = c(-20,25))
ts.plot(Y1[60,,12], main = "REGION 60",
        ylim = c(-20,25))
par(mfrow=c(1,1))

# as we can observe looking at the stationary table, most of the subjects have
# a very good level of stationarity and there are no major differences among these
# numbers. Let's try to get some insights on our original query from these numbers:

# we should also mention that we have very few subjects for getting representative 
# results. However from our small sample qhat we get is the following:


#  Now we'll try to fit several models and then decide which one to include in the final report.

# Let's try to fit a poisson model:

df <- cbind(covariates, stationary = as.numeric(c(stationary.table[,2], stationary.table[,4])))

df.complete <- df[-(1:4),] # 18 variables

# since i suspect that a change in age will affect the current diagnosis, 
# i also include the interaction among these two variables.
# In the other cases interactions should not be necessary.
#I also suspect interaction among lifetime diagnosis and current diagnosis.

pois.model <- glm(formula = stationary ~ age + handedness + current_diagnosis + 
                    lifetime_diagnosis + age*current_diagnosis  
                  # + current_diagnosis*lifetime_diagnosis
                  ,data = df.complete, family = poisson(link = "log")) 
summary(pois.model)

print(xtable(pois.model, type = "latex"), file = "latex_elements/poismodel.tex")


#check the model assumptions

# is the Poisson assumption really met?
# The first possibility that comes to mind for assessing this is to make a 
# histogram of the y’s and see if it looks like the right Poisson distribution.

# Y <- df$stationary
# 
# hist(Y)
# 
# length(df$stationary) #22
# 
# (lambda <- mean(df$stationary))
# 
# tbl <- NULL
# for (k in 0:70) {
#   tbl <- rbind(tbl,c(obs=sum(df$stationary==k),
#                        exp=22*exp(-lambda) * lambda^k / factorial(k)))
#   }
# barplot(tbl[,2])
# 
# # maybe rejected maybe not.... we don't have enough power._.
# 
# ggplot(df, aes(x=Y)) + 
#   geom_histogram(bins = 100)
# 
# mean(Y)
# var(Y)
# we should opt for negative binomial
# https://cran.r-project.org/web/packages/GlmSimulatoR/vignettes/count_data_and_overdispersion.html


# negative binomial regression

NB.model <- glm.nb(stationary ~ age + handedness + current_diagnosis + 
                     lifetime_diagnosis + age*current_diagnosis, 
                   data = df.complete, link = "log")
summary(NB.model)
summary(pois.model)



# checking the validity of the statistic ----------------------------------
plot_series <- function(Y){
  as_tibble(Y) %>%
    mutate(Time = 1:n()) %>%
    pivot_longer(-Time, names_to = "Series", values_to = "vals") %>%
    mutate(Series = factor(Series, levels = colnames(Y))) %>%
    ggplot() +
    geom_line(aes(Time, vals), size =.79, ) + 
    #geom_point(aes(x = Time, y=values), size = 1.5, shape=18, colour="darkred")+
    facet_wrap(facets = vars(Series), ncol = 1) + 
    ylab("") +
    theme_bw() 
}

which(as.numeric(n.stationary[,])==70)
#plot of non stationary TS
ts.plot(t(Y1[c(8,43,66),,8]))
#plot of stationary TS
sample(1:70,1)

ts.plot(t(Y1[c(4,24,35),,8]))

# despite being in a little bit different scale, 
# we ca see there are no major visual differences between stationary and 
# non stationary tseries
# sub 9 vs sub 13
set.seed(4)
ts.plot(Y1[1,,8], ylim = c(-10,15), ylab = "BOLD")
for(i in sample(2:70,10)){
  lines(Y1[i,,8], col=i)
}
#sub 13
set.seed(6)
ts.plot(Y1[1,,12], ylim = c(-10,15), ylab = "BOLD")
for(i in sample(2:70,10)){
  lines(Y1[i,,12], col=i)
}

print(xtable(df[c(8,20,21),], type = "latex"), file = "latex_elements/difftraits.tex")




# investigate the best model ----------------------------------------------
# the idea is to compare different types of regressions 
# based on cross validation 
# regression using age and diagnosis only , and in sepa<rated models.

df <- cbind(covariates, stationary = as.numeric(n.stationary$`n.stationary[, -1]`))

df$diagnosis <- c(rep(NA,4), F,T,T,F,F,F,F,T,F,F,F,T,F,T,T,T,F,F)


df.complete <- df[-(1:4),]


write.csv(df, file="data/regressionADFdataset.csv")


#let's start from the simplest case: OLS

fit.lm <- lm(stationary~ age + diagnosis, df)
summary(fit.lm)
autoplot(fit.lm) 

#find optimal lambda for Box-Cox transformation 
bc <- boxcox(df$stationary ~ df$age + df$diagnosis)
(lambda <- bc$x[which.max(bc$y)])


fit.lm.bc <- lm((stationary^2 -1)/2~ age + diagnosis, df)
summary(fit.lm.bc)
autoplot(fit.lm.bc) 

fit.lm.sqrt<- lm(sqrt(stationary)~ age+diagnosis, df)
summary(fit.lm.sqrt)
autoplot(fit.lm.sqrt)

fit.lm.log<- lm(log(stationary)~ age+diagnosis, df)
summary(fit.lm.log)
autoplot(fit.lm.log)

fit.lm.sq<- lm((stationary)^2~ age+diagnosis, df)
summary(fit.lm.sq)
autoplot(fit.lm.sq)

# Running the LOO-CV on straight lm
df.complete<- df[-(1:4),]
n <- 18 # Number of observations
yhat <- numeric(0)
for (i in 1:n){
  fmi <- lm(stationary~age+diagnosis, data=df.complete[-i,]) # Fit without point i
  yhat[i] <- predict(fmi,df.complete[i,]) # Prediction of y_i
}
sqloss.lm <- sqrt(mean((yhat-df.complete$stationary)^2))
l1loss.lm <- mean(abs(yhat-df.complete$stationary))
sqloss.lm
l1loss.lm

# Running the LOO-CV on sqrt
n <- 18 # Number of observations
yhat <- numeric(0)
for (i in 1:n){
  fmi <- lm(sqrt(stationary)~age+diagnosis, data=df.complete[-i,]) # Fit without point i
  yhat[i] <- predict(fmi,df.complete[i,]) # Prediction of y_i
}
sqloss.lm.sqrt <- sqrt(mean((yhat^2-df.complete$stationary)^2))
l1loss.lm.sqrt <- mean(abs(yhat^2-df.complete$stationary))
sqloss.lm.sqrt # performing worse
l1loss.lm.sqrt

# Running the LOO-CV on square
n <- 18 # Number of observations
yhat <- numeric(0)
for (i in 1:n){
  fmi <- lm((stationary)^2~age+diagnosis, data=df.complete[-i,]) # Fit without point i
  yhat[i] <- predict(fmi,df.complete[i,]) # Prediction of y_i
}
sqloss.lm.sq <- sqrt(mean((sqrt(yhat)-df.complete$stationary)^2))
l1loss.lm.sq <- mean(abs(sqrt(yhat)-df.complete$stationary))
sqloss.lm.sq # performing worse
l1loss.lm.sq


# Running the LOO-CV on log
n <- 18 # Number of observations
yhat <- numeric(0)
for (i in 1:n){
  fmi <- lm(log(stationary)~age+diagnosis, data=df.complete[-i,]) # Fit without point i
  yhat[i] <- predict(fmi,df.complete[i,]) # Prediction of y_i
}
sqloss.lm.log <- sqrt(mean((exp(yhat)-df.complete$stationary)^2))
l1loss.lm.log <- mean(abs(exp(yhat)-df.complete$stationary))
sqloss.lm.log # performing worse
l1loss.lm.log

# Running the LOO-CV on ^5
n <- 18 # Number of observations
yhat <- numeric(0)
for (i in 1:n){
  fmi <- lm((stationary)^130~age+diagnosis, data=df.complete[-i,]) # Fit without point i
  yhat[i] <- predict(fmi,df.complete[i,]) # Prediction of y_i
}
sqloss.lm.p5 <- sqrt(mean(((yhat)^(1/100)-df.complete$stationary)^2))
l1loss.lm.p5 <- mean(abs((yhat)^(1/100)-df.complete$stationary))
sqloss.lm.p5 # performing worse
l1loss.lm.p5

# Let's now try the poisson model

fit.pois <- glm(formula = stationary ~ age + diagnosis   
                # + current_diagnosis*lifetime_diagnosis
                ,data = df.complete, family = poisson(link = "log")) 
summary(fit.pois)

autoplot(fit.pois)

yhat <- numeric(0)
for (i in 1:n){
  fmi <- glm(formula = stationary ~ age + diagnosis   
             # + current_diagnosis*lifetime_diagnosis
             ,data = df.complete[-i,], family = poisson(link = "log")) 
  yhat[i] <- predict(fmi,df.complete[i,]) # Prediction of y_i
}
sqloss.pois <- sqrt(mean(((yhat)^(1/100)-df.complete$stationary)^2))
l1loss.pois <- mean(abs((yhat)^(1/100)-df.complete$stationary))
sqloss.pois # performing way worse
l1loss.pois

#########################################################################
# BRAIN PLOT #
load("data/fMRI-ROI-time-series.RData")
purged.Y <- Y[,1:403, -c(1,21),1]
Y <- purged.Y

##########
ggplot() +
  geom_brain(atlas = dk, position = position_brain(hemi ~ side), show.legend = FALSE) 

##################
brain_labels(dk)
unique(dk$data$label)

write_csv(as.data.frame(brain_labels(dk)), file="data/brainregions")

br <- as.vector(read.csv("data/brainregions"))

any(br$brain.labels.dk.[11]==names(Y[,1,1]))

br$brain.labels.dk.[11]<- "lh-unknown"
br$brain.labels.dk.[46]<- "rh-unknown"

dk$data$label

brain_regions(dk)


s<- NULL
test <- NULL
adf<- matrix(NA, nrow = 22, ncol = 70)
for(i in 1: 22){
  for(j in 1:70){
    s <- as.vector(na.omit(Y[j,,i]))
    test <- adf.test(s, alternative = "stationary")
    adf[i,j] <- ifelse(test$p.value <= 0.01,1,0)
  }
  
  cat(round(i/22*100,0))
  cat("%..")
}


###### create a new toy dataframe

adf.1 <- t(as.data.frame(adf[12,])) # not working, we need it of length 90

dk$data$label# length 90

v <- rep(NA, 70)
for( i in 1:70){
  v[i]<-length(which(dk$data$label == brain_labels(dk)[i] ))
}
which(v==2)

# not working because of different names, let's remove the unknown

brain_labels(dk)[11]
brain_labels(dk)[46]

##### let's try with a simpler approach 

names(which(adf[8,]==0))
brain_labels(dk)[which(adf[8,]==0)]

someData.2= tibble(
  label = brain_labels(dk)[which(adf[8,]==0)],
  unit.root = rep(1, 8)
)

ggplot(someData.2) +
  geom_brain(atlas = dk, 
             position = position_brain(hemi ~ side),
             aes(fill = unit.root), color = "black", show.legend = F) +
  scale_fill_viridis_c(option = "inferno", direction = -1) +
  theme_void() #+
#labs(title = "My awesome title", 
#subtitle = "of a brain atlas plot",
#caption = "I'm pretty happy about this!")



ggseg(someData.2, atlas = dk, 
      colour = "black",
      size = .5, 
      position = "stacked",
      #fill = "white",
      mapping = aes(fill = unit.root), show.legend = T) + 
  scale_fill_distiller(palette = "Set3", na.value = "lightgrey") +
  theme(#legend.justification = c(1),
    legend.position = c("bottom"),
    legend.text = element_text(size = 5)) +
  theme_custombrain(
    plot.background = "white",
    text.colour = "black",
    text.size = 12,
    text.family = "mono"
  ) 


#guides(fill = guide_legend(ncol = 3)) 


