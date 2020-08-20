library(readr)
setwd('e:/rmu_study/project/principal_curve')


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
set.seed(1)
plotSample = data[sample(numPoints, 1.2e4, replace = FALSE),]

plot(plotSample, pch = 20, col = "grey", cex = 0.1)
data1 = as.matrix(data)
sample = as.matrix(plotSample)



source("pgh.R")
source("principal_curve_KDE.R")
library(rgl)


kernel_sigma <- 0.01

#===================================
# compute p0
# ==================================
N <- nrow(sample)
dim <- ncol(sample)
library(iterators)
library(parallel)
library(foreach)
library(doParallel)

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
        res[x - (i-1)*chunk.size,] <-  pgh(sample[x,],sample, N, dim, kernel_sigma)$p
      }# shoudn't run in row, the putin must be matrix
      # return local results
      res
    }
)


stopImplicitCluster()
stopCluster(cl)


p0 <- res2.p


# ================================
# remove the noise point
# ==============================

#threshold <- p0[order(p0)== N*0.04]
p0[p0<0.008*mean(p0)]=0
#p0[p0<threshold] = 0

nz_sample <- sample[p0!=0,]
nz_p0 <- p0[p0!=0]

raw_number <- nrow(nz_sample)
new_number <- (raw_number %/% 6)*6
p0 <- nz_p0[1:new_number]
sample <- nz_sample[1:new_number,]



#======================================
# foreach
#======================================
library(foreach)

library(doParallel)
library(iterators)

library(parallel)
# Real physical cores in the computer
cores <- detectCores(logical=F)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)

# split data by ourselves
N <- nrow(sample)
chunk.size <- N/cores

chunk.size <- new_number/cores
system.time(
  res1.p <- foreach(i=1:cores, .combine='rbind') %dopar%
    {  # local data for results
      res <- matrix(0, nrow=chunk.size, ncol=4)
      # out put the new kernel_sigma in the third coordinate
      for(x in ((i-1)*chunk.size+1):(i*chunk.size)) {
        res[x - (i-1)*chunk.size,] <- principal_curve_KDE(sample, sample[x,], kernel_sigma,1,adaptive = T)
      }# shoudn't run in row, the put in must be matrix
      # return local results
      res
    }
)

stopImplicitCluster()
stopCluster(cl)

pc_projection <- res1.p[,1:2]
kernel <- res1.p[,3]
max(kernel)
min(kernel)
mean(kernel)

#========================
# Plot and save  figure
#=======================
library(ggplot2)
setwd('e:/rmu_study/project/principal_curve/figure')
open3d()

#jpeg(file = '3d.jpg')
plot3d(plotSample$dE,plotSample$E, p0)
plot3d(plotSample$dE,plotSample$E, kernel)
rgl.postscript('03d.pdf','pdf')
#dev.off()

# ===========================
# adaptive2
# ===========================

kernel  <- kernel_sigma*(1/(p0+0.1)^(1/6))
max(kernel)
min(kernel)
mean(kernel)

# chunk.size <- N/cores
cores <- detectCores(logical=F)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)
# split data by ourselves
chunk.size <- N/cores
chunk.size <- new_number/cores
system.time(
  res3.p <- foreach(i=1:cores, .combine='rbind') %dopar%
    {  # local data for results
      res <- matrix(0, nrow=chunk.size, ncol=4)
      # out put the new kernel_sigma in the third coordinate
      for(x in ((i-1)*chunk.size+1):(i*chunk.size)) {
        res[x - (i-1)*chunk.size,] <- principal_curve_KDE(sample, sample[x,], kernel_sigma,1,adaptive = T, p0=p0)
      }# shoudn't run in row, the put in must be matrix
      # return local results
      res
    }
)

stopImplicitCluster()
stopCluster(cl)


iter <- res3.p[,4]

pc_projection <- res3.p[,1:2]

# ====================================

points(pc_projection,col=rgb(0,0,1,0.5),cex=0.5)
