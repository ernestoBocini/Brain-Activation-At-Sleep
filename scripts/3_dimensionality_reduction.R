library(tseries)
library(stats)
library(ggfortify) # for the autoplot
library(MTS) # to use function Eccm
library(factoextra) # to plot cool screeplot
library(corrplot)
library(tidyverse)

load("data/fMRI-ROI-time-series.RData")

# Data cleaning -----------------------------------------------------------

Y <- Y[,,,1]
Y <- Y[,,c(-1,-21)] 
Y <- Y[,-404,]


# six-th subject -----------------------------------------------------------
sub6 <- t(Y[,,7]) 
# let's start the dimensionality reduction analysis by looking at subject 6

glimpse(sub6)

ts.plot(sub6[,2], ylim = c(-22,22), ylab = "BOLD", main = "Subject 6")
for(i in 2:70){
  lines(sub6[,i], col=i)
}

# in order to proceed with the reduction we use PCA

adf.test(sub6[1,])
# stationarity not rejected ... look also at the previous analysis 
# to get more insights 

pcaout <- princomp(x = sub6)
summary(pcaout)

fviz_eig(pcaout, addlabels = TRUE, ylim = c(0, 45))

# we can choose 3 as number of PCs, with a percentage of explained 
# variance of 64%

var <- get_pca_var(pcaout)
# correlation circle (autoplot with a circumference)
fviz_pca_var(pcaout, col.var = "black")

# above is also known as variable correlation plots. It shows the 
# relationships between all variables. It can be interpreted as follow:
#   
# - Positively correlated variables are grouped together.
# - Negatively correlated variables are positioned on opposite sides of the 
#   plot origin (opposed quadrants).
# - The distance between variables and the origin measures the quality of 
#   the variables on the factor map. Variables that are away from the origin 
#   are well represented on the factor map.

# The quality of representation of the variables on factor map is called cos2
# (square cosine, squared coordinates) . You can access to the cos2 as follow:

head(var$cos2, 4)

# can visualize the cos2 of variables on all the dimensions using the corrplot package:

corrplot(var$cos2[,1:3], is.corr=FALSE, rect.lwd = 4, tl.cex = .47, number.cex = .1,
         number.font = 1)
# Note that:
# - A high cos2 indicates a good representation of the variable on the principal component. 
# - A low cos2 indicates that the variable is not perfectly represented by the PCs. 

# For a given variable, the sum of the cos2 on all the principal components is equal to one.
# - The cos2 values are used to estimate the quality of the representation
# - The closer a variable is to the circle of correlations, the better its representation on the factor map (and the more important it is to interpret these components)
# - Variables that are closed to the center of the plot are less important for the first components.

# https://biosakshat.github.io/pca.html#correlation-circle


sub3.pca <- pcaout$scores[,1:3]
toymatrix <- as_tibble(sub3.pca) %>%
  mutate(Time = 1:n()) %>%
  pivot_longer(-Time, names_to = "Series", values_to = "vals") 
values <- toymatrix$vals
values[!(toymatrix$Time==188 | 
       toymatrix$Time==133 |
       toymatrix$Time==355 |
      toymatrix$Time==321 |
        toymatrix$Time==128|
        toymatrix$Time==190 |
        toymatrix$Time==356|
        toymatrix$Time==90 |
        toymatrix$Time==320)
       ] <- NA

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


plot_series(sub3.pca)
df1<-data.frame(Time=1:403,BOLD=sub3.pca[,1])
df2<-data.frame(Time=1:403,BOLD=sub3.pca[,2])
df3<-data.frame(Time=1:403,BOLD=sub3.pca[,3])

ggplot(df1,aes(Time,BOLD))+geom_line(aes(color="First Component"), size = .79)+
  geom_line(data=df2,aes(color="Second Component"), size = .79)+
  geom_line(data=df3,aes(color="Third Component"), size = .79)+
  theme_bw()+
  theme(legend.position="none")


# now that we have explored dimensionality reduction for subject 6
# we can proceed in reducing dimensionality for all subjetcs

# PCA for every subject ---------------------------------------------------

Y.pca <- array(dim = c(3,403,22)
               ,dimnames = list(
               comp = c("Comp.1","Comp.2","Comp.3"),
               time = rownames(as.data.frame(Y[1,,1])),
               subject = rownames(as.data.frame(Y[1,1,])) )
               )

for (i in 1:22) {
  X <- t(Y[,,i])
  pcaout <- princomp(x = X)
  Y.pca[,,i] <- t(pcaout$scores[,1:3])
}

# now that we have a new dataset with multivariate time series with 
# order 3 we can proceed with the identification of the VARIMA model to be used

# see the other R file called Model identification

save(Y.pca, plot_series, file = "data/Y_pca.RData")













