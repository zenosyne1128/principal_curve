#-----------------pgh---------------------------#
# Compute the kernel probability(p), the gradient(g), 
# the Hessian(H) and the sigma inverse(SI) on point x.
#-----------------------------------------------#
pgh <- function(x, data, N, dim, kernel_sigma, adaptive = F, p0 = NaN){
  if(as.matrix(is.nan(p0))[1,]) {
    Kvector <- kernel_vec(x, data, N, dim, kernel_sigma)
    Xvector <- matrix(rep(x,N),nrow = dim)-t(data)
    GX <- Xvector*matrix(rep(Kvector,dim), nrow = dim, byrow = T)
    p <- mean(Kvector)  #is p a vector or a number
    
    if (adaptive == T){
    kernel_sigma <- kernel_sigma*(1/(p+0.1)^(1/2))
    Kvector <- kernel_vec(x, data, N, dim, kernel_sigma)
    Xvector <- matrix(rep(x,N),nrow = dim)-t(data)
    GX <- Xvector*matrix(rep(Kvector,dim), nrow = dim, byrow = T)
    p <- mean(Kvector)
    }
    g <- -rowMeans(GX)/kernel_sigma^2
    H <- 1/(N*(kernel_sigma^4)) * (GX%*%t(Xvector))-diag(dim)*p/(kernel_sigma^2)
    SI <- -(1/p)*H+(1/p)^2*g%*%t(g)
    
  }
  else {
    
    #the mothord is to compute p on every data poin x_i,then regard vector p as input, comput a new density by 
    #using the equation: p = 1/n*\sum 1/h_j*k(x-x_j/h_j),for example, h_j=h/sqrt(p_j)
    h  <- kernel_sigma*(1/(p0+0.1)^(1/6))
    w <- p0/sum(p0)
    Kvector <- kernel_vec(x, data, N, dim, kernel_sigma,h = h)
    Xvector <- matrix(rep(x,N),nrow = dim)-t(data)
    GX <- Xvector*matrix(rep(Kvector,dim), nrow = dim, byrow = T)
    p <- as.numeric(Kvector %*% w)
    p2 <- as.numeric(t(w) %*% (Kvector/(h^dim)) )
    #g <-  (diag(h^(-dim)) %*% w ) %*% (-GX) 
    g <- -rowSums( matrix(rep(w/(h^dim),dim), nrow = dim, byrow = T) * GX)
    H <-(GX*matrix(rep(w/(h^4),dim), nrow = dim, byrow = T)) %*% t(Xvector) - diag(dim)*p2
    SI <- -(1/p)*H+(1/p)^2*g%*%t(g)
  }
  return(list(p=p, g=g, H=H, SI=SI, h=kernel_sigma))
}

#-----------------kernel_vec---------------------------
# Compute the "G_{Sigma_{i}}" vector on point x.
#------------------------------------------------------
kernel_vec <- function(x, data, N, dim, kernel_sigma, h = NaN){
  C <- 1/((2*pi)^(dim/2)*kernel_sigma^dim)
  if (as.matrix(is.nan(h))[1,]){
    dx <- (matrix(rep(x,N),nrow = dim)-t(data))/matrix(rep(kernel_sigma,dim*N),nrow = dim, ncol = N)
    
   } 
  else {
    #C <- 1/((2*pi)^(dim/2)*det(diag(h^dim))) #this is a number
    dx <- (matrix(rep(x,N),nrow = dim)-t(data))/matrix(rep(h,dim),nrow = dim, byrow = T)
    
    }
  K <- C*exp(-0.5*colSums(dx^2)) #not sum up all the colomn, but compute the sum of every colomn
  
  return(K)
}
