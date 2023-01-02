library(MTS)
library(tseries)
library(xtable)

load("data/Y_pca.RData")
#if this gives an error, make sure to run script number 3 before

# VARMA model for subject 2
sub3.pca
Eccm(sub3.pca, maxp=3, maxq = 3)
table2 <- Eccm(sub3.pca, maxp=3, maxq = 3)$pEccm

print(xtable(table2, type = "latex"), file = "pEccm.tex")

# VARMA Models for all the subjects ---------------------------------------
models <- as.data.frame(matrix(nrow = 22, ncol = 2))
rownames(models)<- names(Y.pca[1,1,])
colnames(models)<- c("p","q")

for(s in 1:22){
  eccm <- Eccm(t(Y.pca[,,1]), maxp = 3, maxq = 3)$pEccm
  for( i in 1: nrow(eccm)){
    for(j in 1: ncol(eccm)){
      if(eccm[i,j] >= .02){
        p <- i-1
        q <- j-1
        break
      }
    }
    if(eccm[i,j] >= .02){
      break
    }
    if(i == nrow(eccm) & j==ncol(eccm) & eccm[i,j]<.02){
      p <- "error"
      q <- "error"
    }
    
  }
  models[s,]<- c(p,q)
}


save(models, Y.pca, plot_series, file ="data/models.RData")

# now that I've fitted all the model I can proceed with the
# dimensionality reduction. See next file:
# Outlier detection





