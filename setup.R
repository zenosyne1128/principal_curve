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

#=============================================
source("pgh.R")
source("principal_curve_KDE.R")
library(rgl)
kernel_sigma <- 0.01


pc_projection <- principal_curve_KDE(sample, sample, kernel_sigma,1,adaptive = F)

 points(pc_projection,col=rgb(0,0,1,0.5),cex=0.8)
