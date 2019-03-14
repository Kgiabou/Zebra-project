setwd("C:/Users/Ksq450/Dropbox/Population genetics/Quagga project/DATA") ## setup your worikong directory where you have the data

library(snpMatrix)
library(adegenet)

gen <- read.plink("Quagga_maf0.01") ## read in snpMatrix format
sampl <- read.table("sample info.csv", h=T, sep="\t")
gen2 <- matrix(as.integer(gen@.Data),nrow=nrow(gen@.Data)) ## convert snpMatrix to plain matrix format

sam2 <- sampl[order(sampl$country),] ## order individuals table in order of subspecies
r.names <- as.numeric(rownames(sam2)); r.names ## get the rownames of the converted table
gen3 <- gen2[c(r.names),]

count <- unique(sam2$country) ## get the names for each subspecies
for (i in seq_along(count)) ## count the number of individuals per country
{
newsp <- sam2[sam2$country==count[i],]
print(paste(count[i], nrow(newsp), sep="_"))
}

gen4 <- t(gen3) ## transpose the matrix

table(gen4)

      # 0       1       2       3 
 # 420024   78053  330616 2428659 


# pc <- prcomp(gen4, scale=TRUE, center=TRUE)
# summary(pc) ### Summary of the Principal Components

						# PC1     PC2     PC3     
# Standard deviation     2.29740 1.65425 1.32455 
# Proportion of Variance 0.08378 0.04344 0.02785 
# Cumulative Proportion  0.08378 0.12722 0.15506

names2 <- vector()
count <- unique(sam2$country) 
for (i in seq_along(sam2$country))
{
 ind <- sam2[i,]
 ind2 <- as.character(paste(sam2$country[i], sam2$sample[i], sep="_"))
 names2 <- c(names2, ind2)
}

## Vector that will give a shape for plotting based on the subspecies names ##
Shape <- as.factor(sam2$nominal.subspecies) ## factor of subspecies names
fil <- cbind(species=levels(Shape), c(4,15:20)) ## give a shape number (R reads shapes for plotting in numbers) for each level of the factor
for (k in seq_along(sam2[,5])){ 
  for (j in seq_along(fil[,1])){
    if (sam2[k,5] == fil[j,1]) {
      sam2$shape[k] <- as.numeric(fil[j,2])
    }
  }
}


#### With SNPRelate R Package###
library(SNPRelate)
genofile <- snpgdsBED2GDS("Quagga_maf0.01.bed", "Quagga_maf0.01.fam", "Quagga_maf0.01.bim", "All_sp.gds") ### read plink files and convert to gds file called All_sp.gds
genofile2 <- snpgdsOpen("All_sp.gds") ## read the gds file
sampl2 <- read.table("sample info.csv", h=T, sep="\t")
pcc <-  snpgdsPCA(genofile2, autosome.only=FALSE) ## do a PCA on the gds file
head(round(pcc$varprop*100, 2)) ### Variance importance of the Principal Components


names2 <- vector()
count <- unique(sampl2$country) 
for (i in seq_along(sampl2$country))
 ind <- sampl2[i,]
 ind2 <- as.character(paste(sampl2$country[i], sampl2$sample[i], sep="_"))
 names2 <- c(names2, ind2)
}

## construct a data frame with the pca values of the first two principal components and the names of the samples ##
tab <- data.frame(sample.id = names2,
EV1 = pcc$eigenvect[,1], # the first eigenvector
EV2 = pcc$eigenvect[,2], # the second eigenvector
stringsAsFactors = FALSE)

### Give a shape to each sample based on the subspecies name ##
Shape <- as.factor(sampl2$nominal.subspecies)
fil <- cbind(species=levels(Shape), c(4,15:20))
for (k in seq_along(sampl2[,5])){
  for (j in seq_along(fil[,1])){
    if (sampl2[k,5] == fil[j,1]) {
      sampl2$shape[k] <- as.numeric(fil[j,2])
    }
  }
}

### Give a color to each sample based on the country the sample is from ##
cols <- as.factor(sampl2$country)
fil <- cbind(country=levels(cols), c("skyblue", "purple", "goldenrod", "orange", "grey", "Dark red", "forestgreen")) ## setup color names
for (k in seq_along(sampl2[,3])){
  for (j in seq_along(fil[,1])){
    if (sampl2[k,3] == fil[j,1]) {
      sampl2$cols[k] <- as.character(fil[j,2])
    }
  }
}
## PLot in pdf the results of the PCA ##
pdf("All Species PCA SNPRelate.pdf")
plot(tab$EV2, tab$EV1, xlab="PC2 3.91% ", ylab="PC1 17.21%", col=sampl2$cols, pch=sampl2$shape)
legend("topleft", legend=levels(sampl2$country), pch=16, col=c("skyblue", "purple", "goldenrod", "orange", "grey", "Dark red", "forestgreen"), cex=0.8)
legend("topright", legend=levels(sampl2$nominal.subspecies), pch=c(4,15:20), cex=0.8)
text(0.0, 0.05, "Plains Zebra", cex=0.8)
text(0.0, 0.52, "Grevys zebra", cex=0.8)
dev.off()
