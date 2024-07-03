#TUTORIAL FOR NEMO
##Install NEMO from github
library(devtools)
library(SNFtool) #requires SNF to work
devtools::install_github('Shamir-Lab/NEMO/NEMO')
library(NEMO)

#Load data downloaded from github
#Columns correspond to samples (100 samples/omic) and rows to the features
omic1 = read.table(file='sample_omic1.txt')
omic2 = read.table(file='sample_omic2.txt')

#Create a list of the omics
omic_list <- list(omic1, omic2)

#Affinity graph: first computes the similarity matrix as SNF did, then by using knn 
#it computes a graph, after this it counts the number of times patients share data 
#finally, it divides the affinity matrix by these counts.
affinity_graph=nemo.affinity.graph(omic_list, k=20)

#Estimate number of clusters
num_clusters=nemo.num.clusters(affinity_graph)

#Perform nemo clustering, it does spectral clustering
clustering_nemo=nemo.clustering(omic_list, num.clusters = 3, num.neighbors = 5)

#We can also perform spectral clustering
clustering_spectral = spectralClustering(affinity_graph, num_clusters)

#Use functions from SNF tool to show the different clusters
displayClusters(affinity_graph, clustering_nemo)
displayClustersWithHeatmap(affinity_graph, clustering_nemo)


