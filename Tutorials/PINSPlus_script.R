#TUTORIAL FOR PINSPLUS
#Install from CRAN
library(impute) #Necessary to make the package work
install.packages("PINSPlus")
library(PINSPlus)

#load sample data: rows are the samples and columns are the features
data(AML2004)
data <- as.matrix(AML2004$Gene)

#Perform clustering
#by default it uses kmeans algorithm with kmax=5
#we can also use "pam" and "hclust"
#Perturbation clustering is to cluster a single data type
result <- PerturbationClustering(data=data)
result$k
result$cluster

#Compare results with known subtypes
condition <- seq(unique(AML2004$Group[,2]))
names(condition)=unique(AML2004$Group[,2])
plot(prcomp(AML2004$Gene)$x, col = result$cluster, 
     pch = condition[AML2004$Group[, 2]], main = "AML2004")
legend("topright", legend = paste("Cluster ", sort(unique(result$cluster)), sep = ""),
       fill = sort(unique(result$cluster)))
legend("topleft", legend = names(condition), pch = condition)

#Try with other data
data("KIRC")
dataList <- list(as.matrix(KIRC$GE), as.matrix(KIRC$ME),as.matrix(KIRC$MI))
names(dataList) <- c("GE","ME","MI")
#For multi-omics data we perform another type of clustering
result_subtype <- SubtypingOmicsData(dataList = dataList, clusteringMethod = "kmeans")
result_subtype$cluster1 
result_subtype$cluster2
# Access the cluster assignments
cluster1 <- result_subtype$cluster1
cluster2 <- result_subtype$cluster2

# Plot the clusters for Stage I using a scatter plot
plot(cluster1, pch = 16, col = "blue", main = "Cluster Plot (Stage I)", xlab = "Sample", ylab = "Cluster")
# Plot the clusters for Stage II using a scatter plot
plot(cluster2, pch = 16, col = "red", main = "Cluster Plot (Stage II)", xlab = "Sample", ylab = "Cluster")

#Get clustering outputs
cluster1=result_subtype$cluster1
cluster2=result_subtype$cluster2
a <- intersect(unique(cluster2), unique(cluster1))
names(a) <- intersect(unique(cluster2), unique(cluster1))
a[setdiff(unique(cluster2), unique(cluster1))] <- 
  seq(setdiff(unique(cluster2), unique(cluster1))) + max(cluster1)
colors <- a[levels(factor(cluster2))]

#do survival analysis
coxFit <- coxph(
  Surv(time = Survival, event = Death) ~ as.factor(cluster2),
  data = KIRC$survival,
  ties = "exact"
)
mfit <- survfit(Surv(Survival, Death == 1) ~ as.factor(cluster2), data = KIRC$survival)
plot(
  mfit, col = colors, main = "Survival curves for KIRC, level 2",
  xlab = "Days", ylab = "Survival",lwd = 2
)
legend("bottomright", 
       legend = paste(
         "Cox p-value:", round(summary(coxFit)$sctest[3], digits = 5), sep = ""
       )
)
legend(
  "bottomleft",
  fill = colors,
  legend = paste("Group ", levels(factor(cluster2)), ": ",
                 table(cluster2)[levels(factor(cluster2))], sep =""
  )
)

