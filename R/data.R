#reading data i think
#dont understand this step completely

#install.packages("sparcl")

#library(sparcl)

#read the input data from a CSV file
#x<- read.csv(file = "input.csv", head = FALSE, sep = ",")

#vonvert the data to a matrix
#x <- as.matrix(x)

#hierarchical sparse clustering with permutations
#perm.out <- HierarchicalSparseCluster.permute(x, wbounds = c(1.5, 2:6), nperms = 5)

#create hierarchical sparse cluster using the best w found
#sparsehc <- HierarchicalSparseCluster(dists = perm.out$dists, wbound = perm.out$bestw, method = "complete")

#plot the hierarchical sparse cluster
#plot(sparsehc$hc)

# etach the "sparcl" package
#detach("package:sparcl")
