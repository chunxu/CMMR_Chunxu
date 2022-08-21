#rbiom was developed by Dan at CMMR
#install.packages("rbiom")
library(rbiom)

#load file
biom <- rbiom::read.biom("C:\gao\BCM2021\MILES\Pool912\p912-Goodarzi.metaphlan2.Bacteria.EstCount.txt") 

#normalize data by rarefaction
biom <- rbiom::rarefy(biom)


#rbiom::depth(biom)
rbiom::write.biom(biom, "p912-Goodarzi.metaphlan2.Bacteria.EstCount.rare.txt", format = "tab")
#getwd() #where the file is.
