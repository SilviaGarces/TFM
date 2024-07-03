#TUTORIAL FOR ICLUSTERPLUS AND ICLUSTERBAYES
##Install from bioconductor ##TAKES HOURS TO WORK
library(iClusterPlus)
data(gbm)

#Data pre-processing
##mutation: rows are samples and columns are genes
mut_rate <- apply(gbm.mut,2, mean)
gbm.mut2 <- gbm.mut[,which(mut_rate>0.02)]
gbm.mut2[1:10,1:8]

##Gene expression: use top variable genes
gbm.exp[1:3,1:8]

##Copy number data
gbm.seg[1:3,]

library(GenomicRanges) #Install from bioconductor
data("variation.hg18.v10.nov.2010")
gbm.cn <- CNregions(seg=gbm.seg,epsilon=0,adaptive=FALSE, rmCNV = TRUE, 
                    cnv=variation.hg18.v10.nov.2010[,3:5], frac.overlap=0.5, 
                    rmSmallseg=TRUE, nProbes=5)
dim(gbm.cn)
gbm.cn[1:3,1:5]
gbm.cn <- gbm.cn[order(rownames(gbm.cn)),]
#Check all three datasets are in the same order
all(rownames(gbm.cn)==rownames(gbm.exp))
all(rownames(gbm.cn)==rownames(gbm.mut2))

#Perform clustering- takes a while
fit.single <- iClusterPlus(dt1=gbm.mut2, dt2=gbm.cn, dt3=gbm.exp, type=c("binomial","gaussian","gaussian"),
                           lambda=c(0.04,0.61,0.90), K=2, maxiter=10)

##Tune model
set.seed(123)
#In Windows only 1 CPU - not possible to do it because it takes more than 10 hours if only 1 CPU
for (k in 1:3){
  cv.fit<- tune.iClusterPlus(cpus = 1, dt1 = gbm.mut2, dt2=gbm.cn, dt3=gbm.exp,
                             type=c("binomial","gaussian","gaussian"), K=k, n.lambda = 35,
                             scale.lambda = c(1,1,1),maxiter = 20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

##DO WITH iClusterBayes: it is not faster
fit_bayes <- tune.iClusterBayes(cpus=1,dt1=gbm.mut2, dt2=gbm.cn, dt3=gbm.exp, type=c("binomial","gaussian","gaussian"),
                           K=1:3, n.burnin = 18000, n.draw = 12000, prior.gamma = c(0.5,0.5,0.5), sdev = 0.005, thin = 3)

#Solution would be to do multiple models considering different clusters and check BIC score - the lower the better