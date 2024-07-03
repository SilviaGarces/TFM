#TUTORIAL FOR MIXOMICS
##Install from Bioconductor
library(mixOmics)

#1. UNSUPERVISED MULTIBLOCK (s)PLS
data("breast.TCGA")

##We use 3 types of data, miRNA, mRNA and protein
#Rows are samples and columns are the data type, which is continuous
X1<- breast.TCGA$data.train$mirna
X2<- breast.TCGA$data.train$mrna
X<- list(mirna=X1, mrna=X2)
Y<- breast.TCGA$data.train$protein

##Perform multiblock PLS
block.pls.result <- block.pls(X,Y, design="full")
plotIndiv(block.pls.result)
plotVar(block.pls.result, legend = TRUE)

##Perform multiblock sPLS: we select the number of features to use
list.keepX <- list(mrna=rep(5,2), mirna=rep(5,2))
list.keepY <- c(rep(10,2))

block.spls.result <- block.spls(X,Y, design="full",keepX = list.keepX, keepY = list.keepY)

###Make plots
plotLoadings(block.spls.result, ncomp=1)
plotIndiv(block.spls.result)
plotVar(block.spls.result, title= "PLS Correlations - no cutoff",legend = TRUE)
plotVar(block.spls.result, title= "PLS Correlations - 0.5 cutoff",cutoff=0.5, legend = TRUE)

#2. SUPERVISED - DIABLO:
X1_d <- breast.TCGA$data.train$mirna
X2_d <- breast.TCGA$data.train$mrna
X3_d <- breast.TCGA$data.train$protein
X_d <- list(mirna=X1_d, mrna=X2_d, protein =X3_d)
Y_d <- breast.TCGA$data.train$subtype

##Perform multiblock PLS-DA
result.diablo.tcga <- block.plsda(X_d,Y_d)

##pairwise analysis
list.keepX <- c(25,25)
list.keepY<- c(25,25)
pls1 <- spls(X_d[["mirna"]], X_d[["mrna"]],keepX = list.keepX, keepY = list.keepY)
pls2 <- spls(X_d[["mirna"]], X_d[["protein"]],keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(X_d[["protein"]], X_d[["mrna"]],keepX = list.keepX, keepY = list.keepY)

##plots
plotVar(pls1, cutoff = 0.5,title = "(a) mirna vs mrna",
        legend = c("mirna","mrna"),
        var.names = FALSE, style = 'graphics',
        pch = c(16,17), cex = c(2,2),
        col = c('darkorchid','lightgreen'))
plotVar(pls2, cutoff = 0.5,title = "(a) mirna vs protein",
        legend = c("mirna","protein"),
        var.names = FALSE, style = 'graphics',
        pch = c(16,17), cex = c(2,2),
        col = c('darkorchid','lightgreen'))
plotVar(pls2, cutoff = 0.5,title = "(a) protein vs mrna",
        legend = c("protein","mrna"),
        var.names = FALSE, style = 'graphics',
        pch = c(16,17), cex = c(2,2),
        col = c('darkorchid','lightgreen'))

##correlation
cor(pls1$variates$X, pls1$variates$Y)
cor(pls2$variates$X, pls2$variates$Y)
cor(pls3$variates$X, pls3$variates$Y)

##create design
design <- matrix(0.1, ncol = length(X_d), nrow = length(X_d),
                 dimnames = list(names(X_d), names(X_d)))
diag(design)=0

##create model
basic.diablo.model <- block.splsda(X=X_d, Y=Y_d, ncomp=5, design = design)
###plots
plotIndiv(result.diablo.tcga, legend = TRUE)
plotVar(result.diablo.tcga, legend = TRUE)

##Perform multiblock sPLS-DA
list.keepX_d <- list(mirna = c(16, 17), mrna = c(18,5), protein = c(5, 5))
result.sparse.diablo.tcga <- block.splsda(X_d,Y_d, keepX = list.keepX_d)

##tunning
perf.diablo <- perf(basic.diablo.model, validation = "Mfold",
                    folds = 10, nrepeat = 10)
plot(perf.diablo)

##set ncomp by checking the best one
ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
perf.diablo$choice.ncomp$WeightedVote

##Get variables for final model
test.keepX <- list(mirna= c(5:9, seq(10,18,2), seq(20,30,5)),
                   mrna = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   protein = c(5:9, seq(10, 18, 2), seq(20,30,5)))

#Takes a while - 30 mins aprox
tune.TCGA <- tune.block.splsda(X= X_d, Y= Y_d, ncomp = ncomp, 
                               test.keepX = test.keepX, design = design,
                               validation ='Mfold', folds=5, nrepeat=1,
                               dist = "centroids.dist")
list.keepX <- tune.TCGA$choice.keepX

list.keepX

#Final model using the different tunned parameters
final.diablo.model <- block.splsda(X=X_d, Y=Y_d, ncomp=ncomp,
                                   keepX = list.keepX, design = design)

final.diablo.model$design

#Select features
selectVar(final.diablo.model,block='mrna', comp=1)$mrna$name

###plots
plotDiablo(final.diablo.model,ncomp = 1)
plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE,
          title = 'DIABLO Sample Plots')
plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE,
          title= 'DIABLO')
plotVar(final.diablo.model, var.names = FALSE,
        style = 'graphics', legend = TRUE,
        pch = c(16,17,15), cex = c(2,2,2),
        col = c('darkorchid','brown1','lightgreen'))

circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks = c("chocolate3","grey20", size.labels = 1.5))

#Following plot requires lot of window space to work
network(final.diablo.model, blocks = c(1,2,3),
                      cutoff = 0.4)

plotLoadings(final.diablo.model, ncomp=2, contrib = 'max', method = 'median', trim= TRUE)

#To plot heatmap we can do X11() before
cimDiablo(final.diablo.model) 

#dev.off() #Using this when error in graphics occur
