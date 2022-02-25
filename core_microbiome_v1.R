#this script is used to study core microbiome by extracting clade information from CMMR deliverables and calculate prevalence and average relative abundance
library("ggpubr")
library(caret)
library(dplyr)
library(doParallel)
library(reshape2)
library(ggplot2)
library(data.table) 
library(corrplot)
# Read clade data
my_data <- read.delim2("C:/gao/BCM2018/core_microbiome/Anizome_MARS_1105_1236_1275_Taxa_Summary_Tables/metaphlan3.Bacteria.Phylum.RelAb.txt")
#typeof(my_data) #list
#colnames(my_data) #head line
#nrow(my_data) # number of clade
#ncol(my_data)-2 # number of samples
#length(my_data)

#calculate prevalence and average_relative_abundance for each clade in all samples
res <- matrix (data=NA, nrow=nrow(my_data), ncol=5)
colnames(res) <- c("clade_name","Number_of_samples", "number_of_0", "prevalence", "average_relative_abundance")

for ( i in 1:nrow(res)) {
  res[i,1] <- my_data[i,1]
  res[i,2] <- as.numeric(length(my_data)-2)
  res[i,3] <- as.numeric(length(which(my_data[i,]==0)))
  res[i,4] <- 1-as.numeric(res[i,3])/as.numeric(res[i,2])
  res[i,5] <- 0
  for (j in 2:(ncol(my_data)-1)) {
    res[i,5] <- as.numeric(res[i,5]) + as.numeric(my_data[i,j])
  }
  res[i,5] <- as.numeric(res[i,5])/as.numeric(res[i,2])
} 
res<- as.data.frame(res)
write.csv(res,file = "C:/gao/BCM2018/core_microbiome/phylum-all_samples_result.csv")

################################################################################
#divide samples in groups: by body site (last letter of sample ID)
#colnames(my_data) #head line
#build righ function as in excel
right = function (string, char) {
  substr(string,nchar(string)-(char-1),nchar(string))
}
#test right function
#my_data[1,1]
#right(my_data[2,2],1)

#test columns
#colnames(my_data)[2]
#colnames(my_data)[length(my_data)]
#colnames(my_data)[length(my_data)-1]

#generate list for each body site
my_data_interdigital <- my_data
my_data_groin <- my_data
my_data_dorsal <- my_data
my_data_ear <- my_data
my_data_oral <- my_data
my_data_fecal <- my_data

for (k in 2:(length(my_data)-1)){
  if (right(colnames(my_data)[k],1)!="1")
    my_data_interdigital <- select(my_data_interdigital,-colnames(my_data)[k])
  if (right(colnames(my_data)[k],1)!="2")
    my_data_groin <- select(my_data_groin,-colnames(my_data)[k])
  if (right(colnames(my_data)[k],1)!="3")
    my_data_dorsal <- select(my_data_dorsal,-colnames(my_data)[k])
  if (right(colnames(my_data)[k],1)!="4")
    my_data_ear <- select(my_data_ear,-colnames(my_data)[k])
  if (right(colnames(my_data)[k],1)!="5")
    my_data_oral <- select(my_data_oral,-colnames(my_data)[k])
  if (right(colnames(my_data)[k],1)!="6")
    my_data_fecal <- select(my_data_fecal,-colnames(my_data)[k])
}
#export data by body site
write.csv(my_data_interdigital,file = "C:/gao/BCM2018/core_microbiome/phylum-interdigital.csv")
write.csv(my_data_groin,file = "C:/gao/BCM2018/core_microbiome/phylum-groin.csv")
write.csv(my_data_dorsal,file = "C:/gao/BCM2018/core_microbiome/phylum-dorsal.csv")
write.csv(my_data_ear,file = "C:/gao/BCM2018/core_microbiome/phylum-ear.csv")
write.csv(my_data_oral,file = "C:/gao/BCM2018/core_microbiome/phylum-oral.csv")
write.csv(my_data_fecal,file = "C:/gao/BCM2018/core_microbiome/phylum-fecal.csv")
