library(readr)
setwd('d:/Rfiles/principle_curves/碳的同位素/碳的同位素')


# data = read_csv("./data/raw_data/EventList_58_39.csv")
# data = data[,1:2]
# write_csv(data, path = "./data/experiment1.csv")

data = read_csv("experiment1.csv")

## Apply some transformation
data$dE = (data$dE - min(data$dE)) / (max(data$dE) - min(data$dE))
data$dE = data$dE^(0.8)
data$E = (data$E - min(data$E)) / (max(data$E) - min(data$E))

data$E = data$E^(0.8)

## Sample 50000 points and plot
numPoints = nrow(data) 

plotSample = data[sample(numPoints, 1e4, replace = FALSE),]
plot(plotSample, pch = 20, col = "grey", cex = 0.1)
data1 = as.matrix(data)
sample = as.matrix(plotSample)



source("pgh.R")
source("principal_curve_KDE.R")
library(rgl)

kernel_sigma <- 0.01

#======================================
# foreach
#======================================
library(foreach)
library(doParallel)
# Real physical cores in the computer
cores <- detectCores(logical=F)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
# split data by ourselves
chunk.size <- numPoints/cores


system.time(
  res2.p <- foreach(i=1:cores, .combine='rbind') %dopar%
    {  # local data for results
      res <- matrix(0, nrow=chunk.size, ncol=2)
      for(x in ((i-1)*chunk.size+1):(i*chunk.size)) {
        res[x - (i-1)*chunk.size,] <- principal_curve_KDE(sample, sample[x,], kernel_sigma,1,adaptive = F)
      }# shoudn't run in row, the putin must be matrix
      # return local results
      res
    }
)
stopImplicitCluster()
stopCluster(cl)

#===================================
# adaptive2
# ==================================
N <- nrow(sample)
dim <- ncol(sample)

library(foreach)
library(doParallel)
library(iterators)
library(parallel)
# Real physical cores in the computer 
cores <- detectCores(logical=F)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
# split data by ourselves
chunk.size <- N/cores
system.time(
  res2.p <- foreach(i=1:cores, .combine='rbind') %dopar%
    {  # local data for results
      res <- matrix(0, nrow=chunk.size, ncol=1)
      for(x in ((i-1)*chunk.size+1):(i*chunk.size)) {
        res[x - (i-1)*chunk.size,] <-  pgh(sample[x,],sample , N, dim, kernel_sigma)$p
      }# shoudn't run in row, the putin must be matrix
      # return local results
      res
    }
)
stopImplicitCluster()
stopCluster(cl)

p0 <- res2.p
p0[p0<mean(p0)]=0 # remove the noise point

# ====================================

pc_projection <- principal_curve_KDE(sample, sample, kernel_sigma,1,adaptive = T, p0 = p0)



points(pc_projection,col=rgb(0,0,1,0.5),cex=0.8)
