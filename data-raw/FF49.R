## Load packages ---------------------------------------------------------------

library(tidyverse)
devtools::load_all()

## Real data

XSS <- as.matrix(read.table("data-raw/XSP.txt",header=TRUE))
YS <- as.matrix(read.table("data-raw/YSP.txt",header=TRUE))

XS = matrix(0, nrow = nrow(XSS), ncol = 9)
XS[,1]=XSS[,1]
XS[,2]=XSS[,20]
XS[,3]=XSS[,21]
XS[,4]=XSS[,95]
XS[,5]=XSS[,96]
XS[,6]=XSS[,30]
XS[,7]=XSS[,88]
XS[,8]=XSS[,92]
XS[,9]=XSS[,93]

I = as.numeric(ncol(YS)) # number of assets
P = as.numeric(ncol(XS)) # number of factors
T = as.numeric(nrow(YS)) # number of observations

# Estimationg Beta

RY = 60 # Beta: mounth rolling
RE = T-RY # number of remaing observations
RD = RY-1
TT = as.numeric(nrow(YS))

OLS    = as.list(rep(1:RE))
for (t in 1:RE){
  OLS[[t]] = lm(YS[(t):(t+RD),] ~ XS[(t):(t+RD),])
}

Beta  =  array(0,c(I,RE,(P+1)))
for (t in 1:RE){
  Beta[,t,] = OLS[[t]]$coefficients
}

Beta=Beta[,,-1] # Betas of the Factors

X = Beta[,-RE,] # ajusting for lags

Y = YS[(RY+2):TT,] # ajusting for lags
y = t(Y)

result <- ltm_mcmc(X, y, burnin = 100, iter = 100, K = 3)

dim(result$beta[[1]])

# readr::write_rds(result, "data-raw/result.rds", compress = "xz")
