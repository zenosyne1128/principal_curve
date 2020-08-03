library(readr)
library(ggplot2)
setwd('d:/Rfiles/principle_curves/stream_level/stream_level')

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

plotSample = tib[sample(numPoints, 5e4, replace = T, prob = p),]
plot(plotSample$x,plotSample$y, pch = 20, col = "grey", cex = 0.1)

source("pgh.R")
source("principal_curve_KDE.R")
library(rgl)
kernel_sigma <- 10
                        
sample = as.matrix(plotSample)

pc_projection <- principal_curve_KDE(sample, sample, kernel_sigma,1)
points(pc_projection,col=rgb(0,0,1,0.5),cex=0.3)
