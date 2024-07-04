#Now we want to pre-process the data differently and try to cluster by treatment
##Separation by time = 3 datasets

#Load original datasets:
#metabo <- read.csv("Metabolomics.csv", header=TRUE, sep = "\t")
#proteo_48T <- read.csv("Proteomics_48T.csv", header=TRUE, sep = "\t")
#proteo_96T <- read.csv("Proteomics_96T.csv", header=TRUE, sep = "\t")

#We need to make sure the sample names are in the same order
#proteo_48T$Sample <- gsub("^(M|O)-([0-9]+)-(T[0-2])$", "\\1-\\3-\\2", proteo_48T$Sample)

#Check there are no Nan values, if they are, the Nan value will be the mean of the column
#if(sum(is.na(as.matrix(metabo)))) {
#  for(i in 1:ncol(metabo)){
#    metabo[is.na(metabo[,i]), i] <- mean(as.matrix(metabo[,i]), na.rm = TRUE)
#  }
#}
#if(sum(is.na(as.matrix(proteo_48T)))) {
#  for(i in 1:ncol(proteo_48T)){
#    proteo_48T[is.na(proteo_48T[,i]), i] <- mean(as.matrix(proteo_48T[,i]), na.rm = TRUE)
#  }
#}
#if(sum(is.na(as.matrix(proteo_96T)))) {
#  for(i in 1:ncol(proteo_96T)){
#    proteo_96T[is.na(proteo_96T[,i]), i] <- mean(as.matrix(proteo_96T[,i]), na.rm = TRUE)
#  }
#}

#Put rownames as sample names
#rownames(metabo)<- metabo$Sample
#rownames(proteo_48T)<- proteo_48T$Sample
#rownames(proteo_96T)<- proteo_96T$Sample

#Drop the first column 
metabo <- metabo[-c(1)]
proteo_48T <- proteo_48T[-c(1)]
proteo_96T <- proteo_96T[-c(1)]

#There are samples in metabolomic matrix that do not appear in proteomics 
##Find intersection in rownames
common_samples <- intersect(intersect(rownames(metabo), rownames(proteo_48T)), rownames(proteo_96T))
metabo_commom <- metabo[common_samples, ]
proteo_48T_commom <- proteo_48T[common_samples, ]
proteo_96T_commom <- proteo_96T[common_samples, ]

#Once we have the pre-processing done try and separate by time
metabo_T0 <- subset(metabo_commom, grepl("-T0-",rownames(metabo_commom)))
metabo_T1 <- subset(metabo_commom, grepl("-T1-",rownames(metabo_commom)))
metabo_T2 <- subset(metabo_commom, grepl("-T2-",rownames(metabo_commom)))
prot48T_T0 <- subset(proteo_48T_commom, grepl("-T0-",rownames(proteo_48T_commom)))
prot48T_T1 <- subset(proteo_48T_commom, grepl("-T1-",rownames(proteo_48T_commom)))
prot48T_T2 <- subset(proteo_48T_commom, grepl("-T2-",rownames(proteo_48T_commom)))
prot96T_T0 <- subset(proteo_96T_commom, grepl("-T0-",rownames(proteo_96T_commom)))
prot96T_T1 <- subset(proteo_96T_commom, grepl("-T1-",rownames(proteo_96T_commom)))
prot96T_T2 <- subset(proteo_96T_commom, grepl("-T2-",rownames(proteo_96T_commom)))

#Make list of datasets per time and normalize the data using scale function
omics_T0 <- list(scale(metabo_T0), scale(prot48T_T0),scale(prot96T_T0))
omics_T1 <- list(scale(metabo_T1), scale(prot48T_T1),scale(prot96T_T1))
omics_T2 <- list(scale(metabo_T2), scale(prot48T_T2),scale(prot96T_T2))

# Add names to each component of the lists
names(omics_T0) <- c("Metabo_T0", "Prot48T_T0", "Prot96T_T0")
names(omics_T1) <- c("Metabo_T1", "Prot48T_T1", "Prot96T_T1")
names(omics_T2) <- c("Metabo_T2", "Prot48T_T2", "Prot96T_T2")

#Once we have the datasets separated by time we can proceed with working with the methods and see 
#if they are able to classify by treatment

#Now we want to create a list of the correct subgroups so we can compare the results obtained:
#Extract unique names from all omics
sample_names_T0 <- unique(unlist(lapply(omics_T0, rownames)))
sample_names_T1 <- unique(unlist(lapply(omics_T1, rownames)))
sample_names_T2 <- unique(unlist(lapply(omics_T2, rownames)))

#Create dataframe with names
samples_df_T0 <- data.frame(sample = sample_names_T0, stringsAsFactors = FALSE)
samples_df_T1 <- data.frame(sample = sample_names_T1, stringsAsFactors = FALSE)
samples_df_T2 <- data.frame(sample = sample_names_T2, stringsAsFactors = FALSE)

#Add column to put the treatments
samples_df_T0$treatment <- ifelse(grepl("M-", samples_df_T0$sample), "Mepolizumab",
                          ifelse(grepl("O-", samples_df_T0$sample), "Omalizumab", NA))
#Add column to put the treatments
samples_df_T1$treatment <- ifelse(grepl("M-", samples_df_T1$sample), "Mepolizumab",
                                  ifelse(grepl("O-", samples_df_T1$sample), "Omalizumab", NA))
#Add column to put the treatments
samples_df_T2$treatment <- ifelse(grepl("M-", samples_df_T2$sample), "Mepolizumab",
                                  ifelse(grepl("O-", samples_df_T2$sample), "Omalizumab", NA))

