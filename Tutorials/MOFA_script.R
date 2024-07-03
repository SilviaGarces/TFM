#TUTORIAL FOR MOFA
##Install from bioconductor
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

#Check data: features are in rows and samples in columns (total 200 samples)
data("CLL_data")
#Only use omics data
data <- list(mthy=CLL_data$Methylation,mrna=CLL_data$mRNA, mut=CLL_data$Mutations)
lapply(data,dim)
CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")

#Create MOFA object
MOFA_object <- create_mofa(data)

#Plot data overview - grey indicates missing data: MOFA can handle missing data
plot_data_overview(MOFA_object)

#Specify the options
##Default
data_opts <- get_default_data_options(MOFA_object)

model_opts <- get_default_model_options(MOFA_object)
#1st model
#model_opts$num_factors <- 15
#2nd model
model_opts$num_factors <- 10

train_opts <- get_default_training_options(MOFA_object)
train_opts$seed <- 42

#Start training the model
MOFA_object <- prepare_mofa(MOFA_object, data_options = data_opts,
                            model_options = model_opts, training_options = train_opts)

outfile = file.path(getwd(),"model_2.hdf5") #redirect model to output file
MOFA_object_trained <- run_mofa(MOFA_object, outfile, use_basilisk=TRUE) 
#as it uses python we need to specify use_basilisk=TRUE, if not it will not work

#Overview of the model
slotNames(MOFA_object_trained)
names(MOFA_object_trained@data)
names(MOFA_object_trained@expectations)

##Dimensionality of the factor matrix:
dim(MOFA_object_trained@expectations$Z$group1)

# Add sample metadata to the model
samples_metadata(MOFA_object_trained) <- CLL_metadata

#Generate plots
plot_factor_cor(MOFA_object_trained)
##max_r2= number of factors
plot_variance_explained(MOFA_object_trained, max_r2 = 10)
###Factor 1 is the one with more explained variance for all the different variables

##total variance
plot_variance_explained(MOFA_object_trained, plot_total = T)[[2]]

#Analysis of association: using the metadata
correlate_factors_with_covariates(MOFA_object_trained, 
                                  covariates = c("Gender","died","age"), 
                                  plot="log_pval"
)
###1st model: Factors 13 and 14 have association with 2/3 covariates (died in both)
#1st model:
#plot_factor(MOFA_object_trained, factors = c(13,14), color_by = c("died"))
#2nd model: Factor 5 has 3/3
plot_factor(MOFA_object_trained, factors = c(5), color_by = c("Factor5"))
###weights - factor has explained variance for all variables
###view corresponds to metadata name
plot_weights(MOFA_object_trained, view = "mut", factor = 5, nfeatures = 10,
             scale = T)
#Check top weights and select the one with weight 1
plot_top_weights(MOFA_object_trained, view = "mut",factor = 5,
                 nfeatures = 10, scale = T)

#color_by the one with highest weights
plot_factor(MOFA_object_trained, factors = 5, color_by = "trisomy12",
            add_violin = TRUE, dodge = TRUE)

plot_weights(MOFA_object_trained, view = "mrna", factor= 5, nfeatures = 10)

plot_data_scatter(MOFA_object_trained, view = "mrna", factor = 5, features = 4,
                  sign = "positive", color_by = "trisomy12") + labs(y="RNA expression")

###plot heatmaps without noise
plot_data_heatmap(MOFA_object_trained, view = "mrna", factor = 5, features = 25,
                  denoise = TRUE, cluster_rows= FALSE, cluster_cols= FALSE, 
                  show_rownames= TRUE, show_colnames=FALSE  ,scale="row")

###1st model: REPEAT FOR FACTOR 3 (explained variance high in mutations)
plot_weights(MOFA_object_trained, view = "Mutations", factor = 3, nfeatures = 10,
             scale = T)
#--> most of the weights are around 0, but IGHV has weight 1, this is the main clinical marker for CLL
plot_top_weights(MOFA_object_trained, view = "Mutations",factor = 3,
                 nfeatures = 10, scale = T)

plot_factor(MOFA_object_trained, factors = 3, color_by = "trisomy12",
            add_violin = TRUE, dodge = TRUE)

plot_weights(MOFA_object_trained, view = "mRNA", factor= 1, nfeatures = 10)


plot_data_scatter(MOFA_object_trained, view = "Drugs", factor = 3, features = 4,
                  sign = "positive", color_by = "trisomy12") + labs(y="Drug response")

###plot heatmaps without noise
plot_data_heatmap(MOFA_object_trained, view = "mRNA", factor = 3, features = 25,
                  denoise = TRUE, cluster_rows= FALSE, cluster_cols= FALSE, 
                  show_rownames= TRUE, show_colnames=FALSE  ,scale="row")

##Combination of factors
p <- plot_factors(MOFA_object_trained, factors = c(1,3), color_by = "IGHV",
                  shape_by = "trisomy12", dot_size = 2.5, show_missing = T)

p <- p + geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)
