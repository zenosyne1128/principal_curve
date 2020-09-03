library(readr)
library(ggplot2)
library(rgl)
setwd('e:/rmu_study/project/principal_curve')

data = read.table('data.txt')

data = as.matrix(data)

# -------------------------------------------
# cut the matrix 1/4
# -------------------------------------------


mat = matrix(data = 0,nrow = nrow(data)*ncol(data),ncol = 3)

for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    mat[(i-1)*ncol(data)+j,]=c(i,j,data[i,j])
  }
}
 library(tidyverse)
# #```{r setup, warning = F}
# library(knitr)
# library(rgl)
# library(R.matlab)
# knit_hooks$set(webgl = hook_webgl)
# #```{r test1
plot3d(mat,col = rgb(0,0,1,0.2))
tib = tibble(x=mat[,1],y=mat[,2],h=mat[,3])
tib$h = -(tib$h - max(tib$h)) / (max(tib$h) - min(tib$h))
tib$h = tib$h^(2)
contour(data,nlevels = 45)

# --------------------------
# Sample 50000 points and plot
# --------------------------
numPoints = nrow(tib) 
p = tib$h/sum(tib$h)

set.seed(1)
plotSample = tib[sample(numPoints, 5e4, replace = T, prob = p),]
plot(plotSample$x,plotSample$y, pch = 20, col = "grey", cex = 0.1)

source("pgh.R")
source("principal_curve_KDE.R")

kernel_sigma <- 10
                        
sample = as.matrix(plotSample)
sample = sample[,-3]

# -------------------------
# foreach
# ------------------------
library(iterators)
library(parallel)
library(foreach)
library(doParallel)

# Real physical cores in the computer
cores <- detectCores(logical=F)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)

# split data by ourselves
N <- nrow(sample)
chunk.size <- N/cores

#chunk.size <- new_number/cores
system.time(
  res1.p <- foreach(i=1:cores, .combine='rbind') %dopar%
    {  # local data for results
      res <- matrix(0, nrow=chunk.size, ncol=4)
      # out put the new kernel_sigma in the third coordinate
      for(x in ((i-1)*chunk.size+1):(i*chunk.size)) {
        res[x - (i-1)*chunk.size,] <- principal_curve_KDE(sample, sample[x,], kernel_sigma,1,adaptive = F)
      }# shoudn't run in row, the put in must be matrix
      # return local results
      res
    }
)

stopImplicitCluster()
stopCluster(cl)
pc_projection <- res1.p[,1:2]

points(pc_projection,col=rgb(0,0,1,0.5),cex=0.3)

# ---------------------------
# continuation
# -----------------------------
newmat <- matrix(0,nrow = nrow(data),ncol = ncol(data))
for (i in 1:nrow(pc_projection)) {
  j <- pc_projection[i,1]
  k <- pc_projection[i,2]
  newmat[j,k] <- 1
}

write.csv(newmat,file = 'E:\\rmu_study\\project\\skel2graph3d-matlab\\stream.csv',row.names = F,col.names = F)

library(R.matlab)
path <- ('e:/rmu_study/project/principal_curve/endpiont.mat')
endpiont <- readMat(path)
endpiont <- endpiont[[1]]
endpiont[endpiont==0]=NA
endpiont <- na.omit(endpiont)
points(endpiont,col=rgb(0,0,1,0.5),cex=1)
