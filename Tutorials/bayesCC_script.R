#TUTORIAL FOR BCC
#install bayesCC from Github
library(devtools)
devtools::install_github('ttriche/bayesCC')
library(bayesCC)
data(BRCADATA) #load data

# Fit the model for K=2, ..., 5 to compare the mean adjusted adherence
# This can take a while (30 minutes to 1 hour for each K)
require(parallel)
Ks <- 2:5 # takes a while
names(Ks) <- paste0("K", Ks)
Results <- mclapply(Ks, function(k) {
  bayesCC(BRCAData, K=k, IndivAlpha=T, maxiter=10000)
})

# compute adherences
meanAlphaStar <- c()
upperAlphaStar <- c()
lowerAlphaStar <- c()

for (k in Ks) {
  i <- names(Ks)[which(Ks == k)]
  AlphaRaw <- (Results[[i]]$AlphaVec[,2000:10000] - 1 / k) / (1 - 1 / k)
  AlphaStar <- apply(AlphaRaw, 2, mean)
  meanAlphaStar[i] <- mean(AlphaStar)
  upperAlphaStar[i] <- quantile(AlphaStar, 0.975)
  lowerAlphaStar[i] <- quantile(AlphaStar, 0.025)
}

# Plot AlphaStar (mean adjusted adherence) values for each K
plotCI(Ks, meanAlphaStar, ui = upperAlphaStar, li = lowerAlphaStar, 
       xlab = "Number of clusters (K)", ylab = "Mean adjusted adherence") 
# Maximized when K=3

# Fit model for K=3 (can take 30 min - 1 hour)
K3 <- Results[[2]] #position 2 of results 

# Get alpha values
K3$Alpha

##MAKE PLOTS
# Make principal component plots of clusters
par(mfrow = c(2, 2))

#it compares source clustering (Lbest) ang general clustering (Cbest): meaning that source clustering
#is the scepcific for each dataset while the general is unique for all datasets
# mRNA expression
PCs <- prcomp(t(X1))
plot(PCs$x[, 1], PCs$x[, 2], cex = 0, xlab = "PC 1", ylab = "PC 2", 
     main = expression(paste("GE  (", alpha, " = 0.91)")))
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 1] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 2] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 3] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[1]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 1] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 2] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 3] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[1]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 1] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 2] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 3] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[1]][, 3] == 1, 2], pch = 8)

# DNA methylation
m <- 2
PCs <- prcomp(t(X2))
plot(PCs$x[, 1], PCs$x[, 2], cex = 0, xlab = "PC 1", ylab = "PC 2", 
     main = expression(paste("ME  (", alpha, " = 0.69)")))
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)

# miRNA expression
PCs <- prcomp(t(X3))
plot(PCs$x[, 1], PCs$x[, 2], cex = 0, xlab = "PC 1", ylab = "PC 2", 
     main = expression(paste("miRNA  (", alpha, " = 0.56)")))
m <- 3
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "red", 
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "red", 
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "blue", 
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)

m <- 4
PCs <- prcomp(t(X4))
plot(PCs$x[, 1], PCs$x[, 2], cex = 0, xlab = "PC 1", ylab = "PC 2", 
     main = expression(paste("RPPA  (", alpha, " = 0.70)")))
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "black", 
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "black",
       PCs$x[K3$Cbest[, 1] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "red",
       PCs$x[K3$Cbest[, 2] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 1] == 1, 2], pch = 16)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 2] == 1, 2], pch = 3)
points(PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 1], col = "blue",
       PCs$x[K3$Cbest[, 3] == 1 & K3$Lbest[[m]][, 3] == 1, 2], pch = 8)
