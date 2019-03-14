
library(snpMatrix)
library(adegenet)
library(ape)

dat <- read.plink("Quagga_maf0.01") ### read plink file from .bed extension as snpMatrix
sampl<- read.table("sample info.csv", h=T, sep="\t") ### read the sample info table - convert the xls file first to csv file and then import it in R
sam2 <- sampl[order(sampl$country),] ## order individuals table in order of country name
r.names <- as.numeric(rownames(sam2)); r.names ## get the rownames of the converted table

dat2 <- dat[c(r.names),] ## order the genotypes matrix in country order to match sample info table


dat3 <- matrix(as.integer(dat2@.Data),nrow=nrow(dat2@.Data)) ## convert the snpMatrix to normal matrix with the same number of rows
dat3[dat3==00] <- NA ## convert zeros to NAs
dat3 <- dat3-1 ## convert values of 1,2,3 to 0,1,2

datgenind2 <- df2genind(dat3, sep="") ## convert matrix to genind class


names2 <- vector() ### empty vector
count <- unique(sam2$country) ### get the unique country names from the data frame
for (i in seq_along(sam2$country)) ## construct the names of individuals based on country of sample and sample id (xls file)
{
 ind <- sam2[i,]
 ind2 <- as.character(paste(sam2$country[i], sam2$sample[i], sep="_"))
 names2 <- c(names2, ind2)
}
fac3 <- as.factor(names2) ### create a factor based on each individual


# names2 <- vector()
# subs <- unique(sam2$nominal.subspecies) 
# for (i in seq_along(sam2$nominal.subspecies)) ## construct the names of individuals per subspcesies names
# {
 # ind <- sam2[i,]
 # ind2 <- as.character(paste(sam2$nominal.subspecies[i], sam2$sample[i], sep="_"))
 # names2 <- c(names2, ind2)
# }
# fac5<- as.factor(names2) ### factor per nominal subspecies for each individual

## count the number of individuals per country
for (i in seq_along(count)) 
{
newsp <- sam2[sam2$country==count[i],] ## subset the database based on country name
print(paste(count[i], nrow(newsp), sep="_")) ## print contry and number of samples per country
}
fac4 <- as.factor(c(rep("Botswana", times=8), rep("Kenya2", times=2), rep("Kenya", times=2), rep("Kenya2", times=1), rep("Kenya", times=2), rep("Namibia", times=5), rep("S.Africa", times=1), rep("Tanzania", times=31), rep("Uganda", times=4),
rep("Zambia", times=7))) ## construct factor based on country name

datgenpop3 <- genind2genpop(datgenind2, fac3) ### convert to population genotype class based on fac3 - each individual one "population"
datgenpop4 <- genind2genpop(datgenind2, fac4) ### convert to population genotype class based on fac4 - each country one "population"
datgenpop5 <- genind2genpop(datgenind2, fac5) ### convert to population genotype class based on fac5 - each subspecies one "population"


di3 <- dist.genpop(datgenpop3, method=1) ## estimate genetic distances between individuals
di4 <- dist.genpop(datgenpop4, method=1) ## estimate genetic distances between country level populations
di5 <- dist.genpop(datgenpop5, method=1) ## estimate genetic distances between subspecies level populations

## Construct neighbor joining tree based on genetic distances###
pdf("Nj tree Individuals.pdf")
tr_ind <- bionj(di5) ## based on di3 - individuals level
plot(tr_ind, type="radial", use.edge.length = TRUE, show.tip.label = TRUE, tip.color=rep(c("skyblue", "purple", "goldenrod", "orange", "grey", "Dark red", "forestgreen"), c(8, 7, 5, 1, 31, 4, 7)),
cex=0.6, srt=30)
legend("topleft", legend=levels(sam2$country), pch=16, col=c("skyblue", "purple", "goldenrod", "orange", "grey", "Dark red", "forestgreen"), cex=0.7)
dev.off()
save.image("Genetic Diff.RData")

pdf("NJ country level tree.pdf")
tr_pop <- bionj(di4) ## based on country level populations
plot(tr_pop, type="phylogram", use.edge.length = TRUE, show.tip.label = TRUE, tip.color=c("skyblue", "purple", "black", "goldenrod", "orange", "grey", "Dark red", "forestgreen"),
cex=0.6, srt=30)
legend("topright", legend=levels(sam2$country), pch=16, col=c("skyblue", "purple", "goldenrod", "orange", "grey", "Dark red", "forestgreen"), cex=0.7)
dev.off()
´
# write.table(d, "Country distances.txt", row.names=T, col.names=T, quote=FALSE, sep="\t")
