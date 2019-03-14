
## plotting localities with Zebras colors same as in the PCA analysis ##

library(ggmap)
library(mapproj)
map <- get_map(location = 'Africa', zoom = 4, maptype = "terrain") ## get google map
ggmap(map)

sam <- read.table("sample info.csv", h=T, sep="\t") it in R
sam_count <- sam[order(sam$country),] ## oder individuals based on Country name
count <- unique(sam_count$country) ### get unique country names


names2 <- vector()
for (i in seq_along(sam_count$country)) ## construct sample names based on country and sample id 
{
 ind <- sam_count[i,]
 ind2 <- as.character(paste(sam_count$nominal.subspecies[i], sam_count$sample[i], sep="_"))
 names2 <- c(names2, ind2)
}

## construct a column of different colors for different countries ##
cols_count <- c("skyblue", "purple", "goldenrod", "orange", "grey", "Dark red", "forestgreen") ## colors of each country
final_col <- cbind(count=levels(count), cols_count)
for (k in seq_along(sam_count[,3])){
  for (j in seq_along(final_col[,1])){
    if (sam_count[k,3] == final_col[j,1]) {
      sam_count$cols[k] <- as.character(final_col[j,2])
    }
  }
}

dev.new()
pdf("Mapping Data.pdf")
Country <- sam_count$country
Shape <- as.factor(sam_count$nominal.subspecies) ### construct a factor of numbers based on the different subspecies names
## Pklot with ggplot package ##
ggmap(map)+ geom_point(aes(x = longitude, y = latitude, shape= Shape, color=Country), data = sam_count,
size = 3.5, show.legend=TRUE) + scale_shape_manual(values=c(4,15:20)) + 
scale_color_manual(values=c("skyblue", "purple", "goldenrod", "orange", "grey39", "Dark red", "forestgreen"))
dev.off()
## Save Workspace ##
save.image("Mapping localities Quagga.RData")



