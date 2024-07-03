#TUTORIAL FOR SNF
##Install SNF from CRAN
install.packages("SNFtool")
library(SNFtool)
#Data has patients as rows and the omic data in the columns
#The columns correspond to mRNA expression genes for data 1 and DNA methylation for data 2
data("Data1")
data("Data2")

#Get the true labels (true clusters)
truelabel = c(matrix(1,100,1),matrix(2,100,1));

#Calculate distance between matrices using the pairwise squared Euclidean distance (we can choose other distance measurements)
#We can standardize the data is the values are continuous
Data1_sd = standardNormalization(Data1)
Data2_sd = standardNormalization(Data2)

#This calculates the Euclidean distance
Dist1 = dist2(as.matrix(Data1_sd),as.matrix(Data1_sd)) 
Dist2 = dist2(as.matrix(Data2_sd),as.matrix(Data2_sd))

#Now we convert the distance matrices to similarity matrices using the following function
#K is the number of neighbors where affinities out of the neighborhood are set to 0.
#sigma is the hyperparameter used to do the affinity calculation
#We choose them empirically
W1 <- affinityMatrix(Dist1, K = 20, sigma = 0.5)
W2 <- affinityMatrix(Dist2, K = 20, sigma = 0.5)

#Cluster each matrix individually
displayClustersWithHeatmap(W1, spectralClustering(W1, K = 3))
displayClustersWithHeatmap(W2, spectralClustering(W2, K = 3))

#Clustering of both matrices joined (we can use more matrices as well)
W = SNF(list(W1,W2), 20, 10) #Join matrices

#We can check the optimal clusters before by using the following function
##Range of values from 2 to 5, we can modify the range according to the study
estimateNumberOfClustersGivenGraph(W, 2:5) #We get 2 clusters as the best option

#Perform spectral clustering
C=2 #number of clusters
group=spectralClustering(W,C) #this function does the clustering
displayClusters(W,group)

#Performance of clustering 
SNFNMI=calNMI(group, truelabel)

#Get the labels
M_label=cbind(group,truelabel)
colnames(M_label)=c("spectralClustering","TrueLabel")

#Select colors for each one
M_label_colors=cbind("spectralClustering"=getColorsForGroups(M_label[,"spectralClustering"],
                      colors=c("blue","green")),
                     "TrueLabel"=getColorsForGroups(M_label[,"TrueLabel"],
                      colors=c("orange","cyan")))
#Display the different heatmaps
displayClustersWithHeatmap(W, group, M_label_colors[,"spectralClustering"])
displayClustersWithHeatmap(W, group, M_label_colors)

#Rank NMI scores
NMI_scores <- rankFeaturesByNMI(list(Data1, Data2), W)

## You can get cluster labels for each data point by spectral clustering
labels = spectralClustering(W, C)
plot(Data1, col=labels, main='Data type 1')
plot(Data2, col=labels, main='Data type 2')
