#TUTORIAL FOR MOCLUSTER
##Install mogsa package
library(mogsa)

#load simulation data
#columns are the cell lines and rows are the genes
data("NCI60_4arrays")
tumorType <- sapply(strsplit(colnames(NCI60_4arrays$agilent), split = "\\."),"[",1)
colcode <- as.factor(tumorType)
levels(colcode) <- c("red", "green", "blue", "cyan", "orange",
                     "gray25", "brown", "gray75", "pink")
colcode <- as.character(colcode)

#Perform clustering algorithm
moa <- mbpca(NCI60_4arrays, ncomp=10, k="all", method = "globalScore", option = "lambda1",
             center=TRUE, scale= FALSE, moa= TRUE, svd.solver = "fast", maxiter = 1000)
#ncomp: indicate the latent variables to be calculated
#k=all: no sparsity is introduced
plot(moa, value="eig", type=2)

#permutation test by doing multi-block PCA
r <- bootMbpca(moa, mc.cores = 1, B=20, replace = FALSE, resample = "sample")

#We get the top 3 latent variables are significant, do analysis again with 3 variables
moa_3var <- mbpca(NCI60_4arrays, ncomp=3, k=0.1, method = "globalScore", option = "lambda1",
             center=TRUE, scale= FALSE, moa= TRUE, svd.solver = "fast", maxiter = 1000)
#k=0.1: we only keep 10% variables with non-zero coefficients in the result
#Get scores and compare
scr <- moaScore(moa)
scrs <-moaScore(moa_3var)
#do correlation
diag(cor(scr[, 1:3], scrs))

#Visualize plot in 3D
layout(matrix(1:2, 1, 2))
plot(scrs[, 1:2], col=colcode, pch=20)
plot(scrs[, 2:3], col=colcode, pch=20)
legend("topright", legend = unique(tumorType), col = unique(colcode), pch = 20, cex = 0.8)

#Do hierarchical cluster and plot
hcl <- hclust(dist(scrs))
cls <- cutree(hcl, k=4) #select number of clusters manually
clsColor <- as.factor(cls)
levels(clsColor) <- c("red","blue","orange", "pink")
clsColor<- as.character((clsColor))
heatmap(t(scrs[hcl$order, ]), ColSideColors = colcode[hcl$order], Rowv = NA, Colv=NA)

#Get genes
genes <- moaCoef(moa_3var)
genes$nonZeroCoef$agilent.V1.neg
