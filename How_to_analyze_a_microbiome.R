
#Load the data 
ibd_taxa <- read.csv('C:/gao/BCM2021/How_to_analyze_a_microbiome/ibd_taxa.csv', row.names=1) # taxon x sample matrix
ibd_lineages <- read.csv('C:/gao/BCM2021/How_to_analyze_a_microbiome/ibd_lineages.csv', row.names=1) # lineage x taxon matrix
ibd_metadata <- read.csv('C:/gao/BCM2021/How_to_analyze_a_microbiome/ibd_metadata.csv', row.names=1) # metadata x sample data matrix

head(ibd_taxa)
head(ibd_lineages)
head(ibd_metadata)

dim(ibd_taxa)
dim(ibd_lineages)
dim(ibd_metadata)

#To load the BIOM file into R, you need first to install the biomformat R package.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomformat")

BiocManager::install("phyloseq")

library(biomformat)
library(phyloseq)
ibd_biom=read_biom("C:/gao/BCM2021/How_to_analyze_a_microbiome/ibd_data.biom") #The biomformat package is only reading and writing BIOM files, it does not support visualization. This is the job of phyloseq. 


#convert the BIOM  object into a phyloseq object.
metadata=sample_metadata(ibd_biom)
lineages=observation_metadata(ibd_biom)
counts=biom_data(ibd_biom)
ps=phyloseq(otu_table(as.matrix(counts),taxa_are_rows = TRUE),tax_table(as.matrix(lineages)),sample_data(metadata))

#ps can now be used for further visualisation and analysis 
ps

# assigning column names
colnames(tax_table(ps))=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
colnames(sample_data(ps))=c("SRA_metagenome_name", "Age", "Diagnosis", "Fecal.Calprotectin", "antibiotic", "immunosuppressant", "mesalamine", "steroids")

#example of generated barplot, not from IBD data
A <- rmultinom(5, size=100, prob=c(0.2, 0.4, 0.05, 0.02, 0.15, 0.13, 0.01, 0.04))
B <- rmultinom(5, size=100, prob=c(0.6, 0.25, 0, 0.04, 0.02, 0.06, 0.02, 0))
counts <- cbind(A, B)
barplot(counts)

#visualize IBD data
barplot(as.matrix(ibd_taxa))
#the sample names are not informative and make it hard to read
#cannot see which bacteria are in the samples


# generate number of colours equal to number of phyla
colours <- rainbow(length(unique(ibd_lineages$Phylum)))
color.indices <- match(ibd_lineages$Phylum, unique(ibd_lineages$Phylum))
# generate a vector with color for each taxon
colvec <- colours[color.indices]

#Question 1 Please fill in the number of phyla present in this data set.
length(unique(ibd_lineages$Phylum))

barplot(as.matrix(ibd_taxa), col=colvec)
#legend('topright', fill=colours,legend=unique(ibd_lineages$Phylum))

par(xpd=TRUE, mfrow = c(2,1), mar=c(1, 1, 1, 1))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topright", fill=colours,legend=unique(ibd_lineages$Phylum), cex=.75)
barplot(as.matrix(ibd_taxa), col=colvec,  xaxt='n', yaxt='n', ylim=c(0, 1))
dev.off() # restore default plotting area

# generate number of colours equal to number of class
colours <- rainbow(length(unique(ibd_lineages$Class)))
color.indices <- match(ibd_lineages$Class, unique(ibd_lineages$Class))
# generate a vector with color for each taxon
colvec <- colours[color.indices]

par(xpd=TRUE, mfrow = c(2,1), mar=c(1, 1, 1, 1))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topright", fill=colours,legend=unique(ibd_lineages$Class), cex=.75)
barplot(as.matrix(ibd_taxa), col=colvec,  xaxt='n', yaxt='n', ylim=c(0, 1))
dev.off() # restore default plotting area
length(unique(ibd_lineages$Class))

library(ggplot2)
library(reshape2)
long_data <- cbind(ibd_taxa, ibd_lineages)
long_data <- reshape2::melt(long_data)

ggplot(data=long_data, aes(x=variable, y=value, fill=Phylum))+
geom_bar(position='stack', stat='identity') +
theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + labs(x='Sample', y='Relative abundance')

#Calculating richness using the simulated data 
A <- rmultinom(5, size=100, prob=c(0.2, 0.4, 0.05, 0.02, 0.15, 0.13, 0.01, 0.04))
B <- rmultinom(5, size=100, prob=c(0.6, 0.25, 0, 0.04, 0.02, 0.06, 0.02, 0))
counts <- cbind(A, B)
toy_richness <- colSums(counts != 0)
toy_richness

#Wilcoxon signed-rank test 
groups <- c(rep('A', 5), rep('B', 5))
wilcox.test(toy_richness~groups)

#The Chao 1 estimator
# calculate number of singletons
singletons <- colSums(counts == 1) # 1 0 0 1 1 0 0 1 0 0
doubletons <- colSums(counts == 2) # 0 1 0 0 1 1 1 1 2 2
rares <- singletons / (2*doubletons) # 0 0 0 0 0.5 0 0 0.5 0 0
rares[doubletons == 0] <- 0

chao1 <- toy_richness + rares

chao1

#The Shannon index using generated data set
sums <- apply(counts, 2, sum) # get the column sums (sample sums)
norm_counts <- counts # generate new data frame for normalized data
for (i in 1:ncol(counts)){ # divide per sample by sample total
  norm_counts [,i] <- counts[,i]/sums[i]
}

shannon_div <- function(vector){
    vector <- vector*log(vector)
    # the vector has NA values if the species proportions are 0
    vectorsum <- sum(na.omit(vector))
    return((-1)*vectorsum)
}

shannon_diversities <- apply(norm_counts, 2, shannon_div)

shannon_diversities

#Wilcoxon signed-rank test for Shannon diversity above
wilcox.test(shannon_diversities~c(rep('A', 5), rep('B', 5)))

#calculate Pielou's evenness index
shannon_diversities / log(toy_richness)
#Wilcoxon signed-rank test for Pielou's evenness above
wilcox.test(shannon_diversities / log(toy_richness)~c(rep('A', 5), rep('B', 5)))

#The Simpson index
vec = c(0.4, 0.2, 0.3, 0.1)
1 / sum(vec*vec)


# Diversity Index from IBD data
ibd_metadata$richness = colSums(ibd_taxa != 0)

ibd_metadata$shannon = apply(ibd_taxa, 2, shannon_div)

ibd_metadata$pielou = ibd_metadata$shannon / log(ibd_metadata$richness)

ibd_metadata$simpson = colSums(ibd_taxa * ibd_taxa)

ibd_metadata$invsimpson = 1 / ibd_metadata$simpson

# visualize data
boxplot(shannon~Diagnosis, data=ibd_metadata, ylab='Shannon Diversity')

boxplot(invsimpson~Diagnosis, data=ibd_metadata, ylab='Inverse Simpson', ylim=c(0,40))

#statistics
kruskal.test(shannon~Diagnosis, data=ibd_metadata)

pairwise.wilcox.test(ibd_metadata$shannon, ibd_metadata$Diagnosis, p.adjust.method="BH")

#statistics
colnames(ibd_metadata) 
kruskal.test(shannon~antibiotic, data=ibd_metadata)
pairwise.wilcox.test(ibd_metadata$shannon, ibd_metadata$antibiotic, p.adjust.method="BH")

kruskal.test(shannon~steroids, data=ibd_metadata)
pairwise.wilcox.test(ibd_metadata$shannon, ibd_metadata$steroids, p.adjust.method="BH")

#Making nicer plots
long_data <- reshape2::melt(ibd_metadata)
long_data = long_data[long_data$variable %in% c('shannon', 'pielou', 'simpson', 'invsimpson'),]

ggplot(long_data, aes(x=Diagnosis, y=value, colour=variable)) + geom_boxplot() +
  facet_grid(long_data$variable, scales='free')

#rarefy count data

rarefaction <- function(count.vector=c(),depth=100){
probab=count.vector/sum(count.vector)
return(rmultinom(1,size=depth,prob=probab))
}
# we rarefy an example count vector
count.vector=rmultinom(n=1,size=150,prob=c(0.5,rep(0.1,5)))
rarefaction(count.vector)


#Quantitative metadata Diversity analysis
plot(shannon~Age, data=ibd_metadata, ylab='Shannon')

#correlation
cor(ibd_metadata$shannon, ibd_metadata$Age, method=c("spearman"))
cor.test(ibd_metadata$shannon, ibd_metadata$Age, method=c("spearman"))

#correlation
cor(ibd_metadata$shannon, ibd_metadata$Fecal.Calprotectin, method=c("spearman"))
cor.test(ibd_metadata$shannon, ibd_metadata$Fecal.Calprotectin, method=c("spearman"))

#Comparing species distributions

#relative abundance of Faecalibacterium prausnitzii between IBD and control samples
fp.index=which(rownames(ibd_taxa)=="Faecalibacterium_prausnitzii")
wilcox.test(t(ibd_taxa[fp.index,ibd_metadata$Diagnosis=="Control"]),t(ibd_taxa[fp.index,(ibd_metadata$Diagnosis=="CD" || ibd_metadata$Diagnosis=="UC")]))

#further test to see which condition the proportion of Faecalibacterium tends to be higher:

col1=rgb(0,1,0,0.5)
col2=rgb(1,0,0,0.5)
# determine maximum density
maxy=max(hist(t(ibd_taxa[fp.index,]),breaks="FD",plot=FALSE)$density)
# draw first histogram for control samples
hist(t(ibd_taxa[fp.index,ibd_metadata$Diagnosis=="Control"]),breaks="FD", xlim=c(0,max(ibd_taxa[fp.index,])), ylim=c(0,maxy), prob=TRUE,col=col1, border=col1,xlab="Abundance", main="Histogram")
# draw second histogram for IBD samples
hist(t(ibd_taxa[fp.index,(ibd_metadata$Diagnosis=="CD" || ibd_metadata$Diagnosis=="UC")]),breaks="FD",prob=TRUE,col=col2, border=col2,add=TRUE)
# add a legend
legend("right",legend=c("Control","IBD"), lty = rep(1,2), col = c(col1,col2), merge = TRUE, bg = "white", text.col="black")

#Calculating the Jaccard distance: the intersection of taxa divided by the union of taxa.
# use Simulated data
A <- rmultinom(5, size=100, prob=c(0.2, 0.4, 0.05, 0.02, 0.15, 0.13, 0.01, 0.04))
B <- rmultinom(5, size=100, prob=c(0.6, 0.25, 0, 0.04, 0.02, 0.06, 0.02, 0))
counts <- cbind(A, B)
groups <- c(rep('A', 5), rep('B', 5))

jaccard_sets <- function(a, b){
  taxa <- as.character(c(1:length(a)))
  aset <- taxa[a != 0]
  bset <- taxa[b != 0]
  return(1 - length(intersect(aset, bset))/length(union(aset, bset)))
}
jaccard_sets(counts[,1], counts[,10])

#Calculating Bray-Curtis dissimilarities: take the absolute differences of species abundances across two samples and divide this by the total species abundances.
samples <- cbind(counts[,1], counts[,10]) #take one sample in each group
1 - (2 * sum(apply(samples,1,min)))/sum(samples)
samples

#UniFrac distance:  takes into account relatedness of species to calculate a distance
# no code included
counts[,] #just to show the samples in previous cell

library(vegan)
bray <- vegdist(t(norm_counts))
pcoa.res.sim <- capscale(t(norm_counts)~1, distance='bray', na.action='na.omit')
eigenvalues <- eigenvals(pcoa.res.sim)

eigenvalues

eigenvalues/sum(eigenvalues) #What percentage of variation is contained in the eigenvector?

plot(pcoa.res.sim)

#remove labels 
plot(pcoa.res.sim$CA$u[,c(1,2)],
     xlab=paste(c('PCoA1', round(eigenvalues[[1]],2))), ylab=paste(c('PCoA2', round(eigenvalues[[2]],2))))

#Please specify the name of the parameter that takes the colvec vector and then colours the samples in the plot function.
plot(pcoa.res.sim$CA$u[,c(1,2)],col = c(rep("#FF0000FF", 5), rep("#00FFFFFF", 5)), xlab=paste(c('PCoA1', round(eigenvalues[[1]],2))), ylab=paste(c('PCoA2', round(eigenvalues[[2]],2))))

#visualize the 2nd axis as the x-axis and the 3rd axis as the y-axis
plot(pcoa.res.sim$CA$u[,c(2,3)],col = c(rep("#FF0000FF", 5), rep("#00FFFFFF", 5)), xlab=paste(c('PCoA1', round(eigenvalues[[1]],2))), ylab=paste(c('PCoA2', round(eigenvalues[[2]],2))))

#Non-metric multidimensional scaling 
nmda.res.sim <- metaMDS(t(norm_counts), distance="bray", k=2, trymax=50)
stressplot(nmda.res.sim)

plot(nmda.res.sim)

colvec = c(rep("#FF0000FF", 5), rep("#00FFFFFF", 5))
ordiplot(nmda.res.sim,type="n")
orditorp(nmda.res.sim, display="sites", col=colvec)
points(nmda.res.sim, display="species", col='green')

'''
To be continued....
'''

