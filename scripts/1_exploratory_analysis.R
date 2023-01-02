#loading useful libraries
library(tseries)

library(xtable)


# load the datasetand the covariates:
load("data/fMRI-ROI-time-series.RData") 
covariates <- read.table("data/subj-covariates.txt", header = T)
covariates$type_current_diagnosis[is.na(covariates$type_current_diagnosis)] <- "-"
covariates$type_lifetime_diagnosis[is.na(covariates$type_lifetime_diagnosis)] <- "-"

#export table for latex
print(xtable(covariates, type = "latex"), file = "latex_elements/covariates.tex")

# explore the data: Dynamic activity time series of each brain region for the second subject.
ts.plot(Y[1,,2,1], ylim = c(-10,15), ylab = "BOLD")
for(i in 2:70){
  lines(Y[i,,2,1], col=i)
}

# analysing the missing values for the scan and the rescan subject-wise

Y1 <- Y[,,,1]
Y2 <- Y[,,,2]

missing.1<- NULL 
missing.rate  <- NULL
missing.2 <- NULL
missing.rate2 <- NULL

for(i in 1:24){
  missing.1[i] <- sum(is.na(Y1[,,i]))
  missing.rate[i]<- paste(round(sum(is.na(Y1[,,i]))/(404*70) *100,4),"%")
  missing.2[i] <- sum(is.na(Y2[,,i]))
  missing.rate2[i]<- paste(round(sum(is.na(Y2[,,i]))/(404*70) *100,4),"%")
}

table.NA <- cbind(subject = colnames(Y1[,1,]), 
                  missing.1, 
                  missing.2,
                  missing.rate.1 = missing.rate,
                  missing.rate.2 = missing.rate2
                  )
rownames(table.NA) <- table.NA[,1]
table.NA <- table.NA[,-1]
#View(table.NA)
# looking at the table we can observe that the number of missing values it's the same for the 
# following subjects:
which(table.NA[,1]=="70")
# let's investigate further this fact:

toymatrix <- Y1[,,c(2,4,5,14,20)]

index.NA <- matrix( nrow = 70, ncol = 5)

for(j in 1:5){
  
  for(i in 1:70){
    if(sum(is.na(toymatrix[i,,j])) != 1){
      index.NA <- "ERROR"
      break
    }else{
      index.NA[i,j] <- as.numeric(which(is.na(toymatrix[i,,j]))) 
    }
  }
  
}
index.NA
# we can immediately see that the missing value is consistently allocated in the last 
# observation. We proceed by removing it.

Y1[3,1,2] <- NA

#export table for latex
xtable(table.NA, type = "latex", file = "latex_elements/tableNA.tex")
print(xtable(table.NA, type = "latex"), file = "latex_elements/latex_elementableNA.tex")

# since the number of NA is sensibly smaller than for the rescan, 
# I decided to drop the second scanning, and I'll proceed with the 
# analysis of the first one. An possible improvement for the global project
# could be to consider the second rescan as well.

# Moreover, the number of missing values is critical only for
# subject_1 (100%) and for subject_21 (100%), while for the others 
# is always smaller than or equal to 0.25%. Once again, we remove these subjects
# which do not bring any information.

###########################################
# Take home message:
# Only consider first scan
# Remove last observation of all subjects
# Remove subject 1 and 21


