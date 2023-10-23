library(WGCNA)
library(cluster)


path <- "./discovery3KAgeGapPerGOSub.csv"
datExpr <- read.table(path, sep=",", header=TRUE)
rownames(datExpr) <- datExpr$X
datExpr$X <-  NULL
datExpr$Age <- NULL

trueModule <- "./analysis/GO_AncestorReTag3000.csv"
trueModule <- read.table(trueModule, sep=",", header=TRUE)
trueModule$GO <- sub("GO:", "GO.", trueModule$GO)
trueModule <- subset(trueModule, GO %in% colnames(datExpr))
trueModuleColor  <- trueModule$Color

# step1 scale free topology
# here we define the adjacency matrix using soft thresholding with beta=6
# beta 2-20
ADJ1=abs(cor(datExpr,use="p"))^6

# When you have relatively few genes (<5000) use the following code
# the sum value of each column
k=as.vector(apply(ADJ1,2,sum, na.rm=T))

# Plot a histogram of k and a scale free topology plot
# Check whether the adjacency is a scale free topology
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

# step2 hierarchical clustering
# Turn adjacency into a measure of dissimilarity
# Adjacency can be used to define a separate measure of similarity, the Topological Overlap Matrix(TOM)
dissADJ=1-ADJ1
dissTOM=TOMdist(ADJ1)
collectGarbage()

# step3 define models
# Module definition using the topological overlap based dissimilarity
# Calculate the dendrogram
hierTOM = hclust(as.dist(dissTOM),method="average");


cutHeightList <-  c(0.99, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999)
minClusterSizeList <- c(10, 20, 30)

# Create a new directory
directory_name <- "discovery3KGOClusterSub"
dir.create(directory_name)
# Check if the directory was created successfully
if (file.info(directory_name)$isdir) {
  cat("Directory created:", directory_name, "\n")
} else {
  cat("Failed to create the directory:", directory_name, "\n")
}
path_out <- "./discovery3KGOClusterSub/"
randIndexResults <- data.frame(CutHeight = numeric(), MinClusterSize = integer(), RandIndex = numeric()) 
for(cutHeight in cutHeightList) {                                            # Head of first for-loop
  for(minClusterSize in minClusterSizeList) {                                          # Head of nested for-loop
        colorDynamicHybridTOM = cutreeDynamic(
            hierTOM, 
            method = "hybrid", 
            distM= dissTOM , 
            cutHeight = cutHeight,
            minClusterSize = minClusterSize, 
            deepSplit=2, 
            pamRespectsDendro = FALSE
        )
        trueModule$Label <- colorDynamicHybridTOM
        tabDynamicHybridTOM =table(colorDynamicHybridTOM,trueModuleColor)
        RandIndex <- randIndex(tabDynamicHybridTOM ,adjust=F)
        randIndexResults <- rbind(randIndexResults, data.frame(cutHeight = cutHeight, minClusterSize = minClusterSize, randIndex = RandIndex))
        write.csv(trueModule, paste(path_out, "trueModuleCluster" ,cutHeight, minClusterSize,  ".csv", sep = "_"), row.names=TRUE)
  }
}
print(randIndexResults)
write.csv(randIndexResults, paste(path_out, "moduleRandIndex.csv", sep = ""), row.names=TRUE) 
