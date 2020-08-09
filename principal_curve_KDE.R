#--------------principal_curve_KDE--------------------------------------#
# Subspace Constrained Mean Shift(SCMS)
# to obtain principal curves/surfaces based on kernel density estimation(KDE);
# X: data points we've got(used for KDE);
# Xinit: points to be projected;
# kernel_sigma: the Gaussian kernel bandwidth;
# targetdim: the dimension of target subspace.
#-----------------------------------------------------------------------#
principal_curve_KDE <- function(X, Xinit, kernel_sigma, targetdim, adaptive = F, p0 = NaN){
  point <- Xinit 
  N <- nrow(X)
  dim <- ncol(X)
  threshold <- 1e-10 #tolerance
  pc <- matrix(0, nrow = N, ncol = dim) 
  
  for (ind in 1:nrow(point)) {
    flag <- 0 
    pghlist <- pgh(point[ind, ], X, N, dim, kernel_sigma, adaptive, p0)
    point_pc1 <- point[ind, ]
    eigenlist <- eigen(pghlist$SI)
    ConstrainedSpace <- eigenlist$vectors[, 1:(dim-targetdim)]#[,dim:(dim-targetdim)]
    ConstrainedSpace <- as.matrix(ConstrainedSpace)
    gra <- pghlist$g
    H <- pghlist$H
    ##########################################################################
    if(abs(t(gra)%*%H%*%gra/(norm(gra,"2")*norm(t(gra)%*%H,"2")))>1-1e-10){
      ##########################################################################
      flag <- 1 #indicates the point already on pricipal curves/surfaces
      pc[ind, ]= point_pc1
    }
    if(!flag){
      for (a in 1:20) {
        G <- kernel_vec(point_pc1, X, N, dim, kernel_sigma, p0)
        num1 <- rowSums(matrix(rep(G, dim), nrow = dim, byrow = T)*t(X))
        den1 <- sum(G)
        #######################################################################
        pghlist <- pgh(point_pc1, X, N, dim, kernel_sigma, adaptive, p0)
        eigenlist <- eigen(pghlist$SI)
        ConstrainedSpace <- eigenlist$vectors[, 1:(dim-targetdim)]
        ConstrainedSpace <- as.matrix(ConstrainedSpace)
        #######################################################################
        point_pc1_old <- point_pc1
        for (c in 1:ncol(ConstrainedSpace)) {
          direction <- ConstrainedSpace[, c]
          point_pc1 <- point_pc1 + direction%*%(t(direction)%*%(num1/den1-point_pc1))
        }
        if(sum(abs(point_pc1_old-point_pc1)<threshold)==dim){
          break
        }
      }
      
      pc[ind, ]= point_pc1
    }
  }
  return(pc)
}