library(marima)
library(forecast)
library(tidyverse)
library(tseries)
library(xtable)
load("data/Y_pca.RData")
load("data/models.RData")

# FUNCTION out.det --------------------------------------------------------
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
indicator <- function(index){
  ar <- seq(1,models$p[index], by = 1)
  ma <- seq(1,models$q[index], by = 1)
  kvar <- 3
  
  indicator<- define.model(kvar = kvar, #how many variables (ts)?
                            ar = ar, 
                            ma=ma, 
                            rem.var = NULL,  # remove some variables?
                            indep = NULL # are variables independent?
  )
  return(indicator)
}

out.det <- function(subj, index = NA, indicator = NULL, kvar = 3, L=5){
  if(is.null(indicator)){
    indicator <- indicator(index)
  }
  
  marima <- quiet(marima(t(subj), means=1,
                   ar.pattern = indicator$ar.pattern, 
                   ma.patter = indicator$ma.pattern,
                   Check = F,
                   #Plot = "log.det", # plot the log of the determinant
                   penalty = 0))

  #### OUTLIER DETECTION ALGORITHM ####
  est.residuals <- marima$residuals
  sigma.res <- marima$resid.cov
  est.ar <- pol.mul(pol.inv(marima$ma.estimates, L=L), marima$ar.estimates, L=L)
  
  weight.A <- as.data.frame(matrix(NA,nrow = nrow(subj), ncol=3))
  weight.I <-  as.data.frame(matrix(NA,nrow = nrow(subj),ncol=3))
  
  J.I <- rep(NA, nrow(subj))
  J.A <- rep(NA, nrow(subj))
  C.I <- rep(NA,nrow(subj))
  C.A <- rep(NA, nrow(subj))
  weight.A.denominator <- -solve(t(est.ar[,,1])%*%solve(sigma.res)%*%(est.ar[,,1]) + 
                                   t(est.ar[,,2])%*%solve(sigma.res)%*%est.ar[,,2]+
                                   t(est.ar[,,3])%*%solve(sigma.res)%*%est.ar[,,3]+
                                   t(est.ar[,,4])%*%solve(sigma.res)%*%est.ar[,,4]+
                                   t(est.ar[,,5])%*%solve(sigma.res)%*%est.ar[,,5]+
                                   t(est.ar[,,6])%*%solve(sigma.res)%*%est.ar[,,6])
  for(h in 7:(nrow(subj)-5)){
    weight.I[h,]<- est.residuals[,h]
    J.I[h] <- t(est.residuals[,h])%*%solve(sigma.res)%*%est.residuals[,h]
    max.I <- as.numeric(which(est.residuals[,h]==max(est.residuals[,h])))
    C.I[h] <- est.residuals[max.I,h]/sqrt(sigma.res[max.I,max.I])
    
    weight.A.numerator <- t(est.ar[,,1])%*%solve(sigma.res)%*%est.residuals[,h]
    for(i in 2:6){
      t = t(est.ar[,,i])%*%solve(sigma.res)%*%est.residuals[,h+i-1]
      weight.A.numerator <- weight.A.numerator + t
    }
    
    weight.A[h,] <- weight.A.denominator%*%weight.A.numerator
    
    sigma.est <- -weight.A.denominator
    
    J.A[h] <- as.numeric(weight.A[h,])%*%solve(sigma.est)%*%t(t(as.numeric(weight.A[h,])))
    max.A <- as.numeric(which(weight.A[h,]==max(weight.A[h,])))
    C.A[h] <- weight.A[h,max.A]/sqrt(sigma.est[max.A,max.A])
  }
  
  J.scores <- cbind(J.I = as.vector(na.omit(J.I)),J.A=as.vector(na.omit(J.A)),
                    C.I = as.vector(na.omit(C.I)), C.A=as.vector(na.omit(C.A)))
  return(list(Scores = J.scores, weight.A = weight.A,weight.I =  weight.I))
}

# SUBJECT 2 -------------------

subj2 <- t(Y.pca[,,1])
plot_series(subj2)
ar <- seq(1,models$p[1], by = 1)
ma <- seq(1,models$q[1], by = 1)
kvar <- 3

indic <- indicator(1)

out.scores.sub2 <- out.det(subj = subj2, index = 1, indicator = indic, kvar = 3, L=5)
out.scores.sub2$Scores
head(out.scores.sub2)
outscores <- head(out.scores.sub2$Scores[order(out.scores.sub2$Scores[,1], decreasing = T),])
outscores <- cbind(order(out.scores.sub2$Scores[,1], decreasing = T)[1:6], outscores)
print(xtable(outscores, type = "latex"), file ="latex_elements/outscores.tex")

# plot of outliers for subject 2--------------------------------------------------------
Scores <- out.scores.sub2$Scores
max(Scores[,1])
which(Scores[,2] == max(Scores[,2]))
Scores[order(Scores[,1], decreasing = T),1]

toymatrix <- as_tibble(subj2) %>%
  mutate(Time = 1:n()) %>%
  pivot_longer(-Time, names_to = "Series", values_to = "vals")
toymatrix$vals2 <- toymatrix$vals

i.o. <- toymatrix$vals[c(
  which(Scores[,1]>=16.19)*3-2 + 6*3,
  which(Scores[,1]>=16.19) *3-1 + 6*3,
  which(Scores[,1]>=16.19) *3 + 6*3)] 
toymatrix$vals <- NA
toymatrix$vals[c(
  which(Scores[,1]>=16.19)*3-2 + 6*3,
  which(Scores[,1]>=16.19) *3-1 + 6*3,
  which(Scores[,1]>=16.19) *3 + 6*3)] <- i.o.

a.o. <- toymatrix$vals2[c(
  which(Scores[,2]>=16.21) *3-2 + 6*3,
  which(Scores[,2]>=16.21) *3-1 + 6*3,
  which(Scores[,2]>=16.21) *3 + 6*3)] 
toymatrix$vals2 <- NA
toymatrix$vals2[c(
  which(Scores[,2]>=16.21)*3-2 + 6*3,
  which(Scores[,2]>=16.21)*3-1 + 6*3,
  which(Scores[,2]>=16.21)*3 + 6*3)] <- a.o.

    as_tibble(subj2) %>%
      mutate(Time = 1:n()) %>%
      pivot_longer(-Time, names_to = "Series", values_to = "vals") %>%
      mutate(Series = factor(Series, levels = colnames(subj2))) %>%
      ggplot() +
      geom_line(aes(Time, vals)) + 
      facet_wrap(facets = vars(Series), ncol = 1) + 
      geom_point(aes(x = Time, y = toymatrix$vals , col="red"))+
      geom_point(aes(x = Time, y = toymatrix$vals2 , col="blue"))+
      ylab("") +
      theme_bw()



# SUBJECT 3 ---------------------------------------------------------------

subj3 <- t(Y.pca[,,2])
ar <- seq(1,models$p[2], by = 1)
ma <- seq(1,models$q[2], by = 1)
plot_series(subj3)

indic <- define.model(kvar = kvar, #how many variables (ts)?
                          ar = ar, 
                          ma=ma, 
                          rem.var = NULL,  # remove some variables?
                          indep = NULL # are variables independent?
)

out.scores.sub3 <- out.det(ar,ma,subj = subj3, indicator=indic)

Scores <- out.scores.sub3$Scores

which(Scores[,1]>=16.19)
which(Scores[,2]>=16.21)
Scores[order(Scores[,1], decreasing = T),1]
which(Scores[,1]== max(Scores[,1]))


# Global Results ----------------------------------------------------------


global_results <- function(subject=1, o=1){
  subj2 <- t(Y.pca[,,subject])
  ar <- seq(1,models$p[subject], by = 1)
  ma <- seq(1,models$q[subject], by = 1)
  kvar <- 3
  
  indicator <- define.model(kvar = kvar, #how many variables (ts)?
                            ar = ar, 
                            ma=ma, 
                            rem.var = NULL,  # remove some variables?
                            indep = NULL # are variables independent?
  )
  
  out.scores.sub2 <- out.det(subj = subj2, indicator = indicator)$Scores
  #out.scores.sub2 <- cbind(rownames(subj2)[7:398], as.numeric(out.scores.sub2))
  head(out.scores.sub2)
  outscores <- head(out.scores.sub2[order(out.scores.sub2[,o], decreasing = T),])
  outscores <- cbind(order(out.scores.sub2[,o], decreasing = T)[1:6] + 6, outscores)
  return(outscores)
}

global_results(subject = 3)

psi.fun <- function(index){
  library(marima)
  subj <- t(Y.pca[,,index])
  indicator <- indicator(index)
  marima <- quiet(marima(t(subj), means=1,
                         ar.pattern = indicator$ar.pattern, 
                         ma.patter = indicator$ma.pattern,
                         Check = F,
                         #Plot = "log.det", # plot the log of the determinant
                         penalty = 0))
  L <- 5
  psi <- pol.mul(pol.inv(marima$ar.estimates, L=L), marima$ma.estimates, L=L)
  return(psi)
}

out_algorithm <- function(subject.number, Y.pca) {
  subject.n <- subject.number - 1
  subj <- t(Y.pca[, , subject.n])
  outliers <- out.det(subj = subj, index = subject.n)
  psi <- psi.fun(subject.n)
  
  out.unit <- matrix(NA, nrow = 20, ncol = 3)
  
  for(iterations in 1:20) {
    gsI <- global_results(subject.n, 1)[1, 1:2]
    gsA <- global_results(subject.n, 2)[1, c(1, 3)]
    gsIc <- global_results(subject.n, 3)[1, c(1, 4)]
    gsAc <- global_results(subject.n, 4)[1, c(1, 5)]
    
    if (gsI[2] > gsA[2] & gsI[2] > 32.24) {
      out.unit[iterations, 1] <- as.integer(gsI[1])
      out.unit[iterations, 2] <- "J_IO"
      out.unit[iterations, 3] <- gsI[2]
      
      for (i in 0:5) {
        Y.pca[, gsI[1] + i, subject.n] <-
          subj[gsI[1] + i,] - t(psi[, , i + 1] %*% as.numeric(outliers$weight.I[gsI[1] + i,]))
      }
      
      cat(paste0(iterations, ".."))
    }
    if (gsI[2] < gsA[2] & gsA[2] > 32.66) {
      out.unit[iterations, 1] <- gsA[1]
      out.unit[iterations, 2] <- "J_AO"
      out.unit[iterations, 3] <- gsA[2]
      
      Y.pca[, gsA[1], subject.n] <-
        subj[gsA[1],] + as.numeric(outliers$weight.A[gsA[1],])
      
      cat(paste0(iterations, ".."))
      
    }
    if (gsI[2] < 32.66 & gsA[2] < 32.66) {
      if (gsIc[2] > 3.90 & gsAc[2] < 4.56) {
        out.unit[iterations, 1] <- gsIc[1]
        out.unit[iterations, 2] <- "C_IO"
        out.unit[iterations, 3] <- gsIc[2]
        for (i in 0:5) {
          Y.pca[, gsIc[1] + i, subject.n] <-
            subj[gsIc[1] + i,] - t(psi[, , i + 1] %*% as.numeric(outliers$weight.I[gsIc[1] +
                                                                                     i,]))
        }
        
        cat(paste0(iterations, ".."))
        
      }
      
      if (gsIc[2] < 3.90 & gsAc[2] > 4.56) {
        out.unit[iterations, 1] <- gsAc[1]
        out.unit[iterations, 2] <- "C_AO"
        out.unit[iterations, 3] <- gsAc[2]
        Y.pca[, gsAc[1], subject.n] <-
          subj[gsAc[1],] + as.numeric(outliers$weight.A[gsAc[1],])
        
        cat(paste0(iterations, ".."))
        
      }
      
      if (gsIc[2] > 3.90 & gsAc[2] > 4.56) {
        out.unit[iterations, 1] <- gsAc[1]
        out.unit[iterations, 2] <- "C_AO"
        out.unit[iterations, 3] <- gsAc[2]
        Y.pca[, gsAc[1], subject.n] <-
          subj[gsAc[1],] + as.numeric(outliers$weight.A[gsAc[1],])
        
        cat(paste0(iterations, ".."))
        
      }
      
      if (gsIc[2] < 3.90 & gsAc[2] < 4.56) {
        break
      }
    }
    
  }
  return(out.unit)
}

save(plot_series, out.det, indicator, quiet, global_results, Y.pca, models, psi.fun, out_algorithm, file = "data/outdet.RData")


load("data/outdet.RData")





