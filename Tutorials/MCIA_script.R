#TUTORIAL FOR MCIA
##Install from bioconductor
library(omicade4)

#Load data from package
data("NCI60_4arrays")
##Data needs to have the same columns for all the datasets but not the same number of rows
#Samples in columns and variables in rows --> same samples in all omics 

#This generates the MCIA model
mcoin <- mcia(NCI60_4arrays)

#Plot results
plot.mcia(mcoin, sample.lab=FALSE, df.col=4:7) 
plotVar(mcoin, var=NA, bg.var.col=1:4, var.lab=TRUE)
plotVar(mcoin, var=c("SPOPL", "CAPN2", "SNX8"), #consider some variables
        df=1:4, var.lab=TRUE, var.col=c("red", "green", "blue"))

#Select cancer type
cancer_type <- colnames(NCI60_4arrays$agilent)
cancer_type <- sapply(strsplit(cancer_type, split="\\."), function(x) x[1])

#Plot model considering the cancer types
plot(mcoin, axes=1:2, phenovec=cancer_type, sample.lab=FALSE, df.color=1:4)

