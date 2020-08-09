#-----------------pgh---------------------------#
# Compute the kernel probability(p), the gradient(g), 
# the Hessian(H) and the sigma inverse(SI) on point x.
#-----------------------------------------------#
pgh <- function(x, data, N, dim, kernel_sigma, adaptive = F, p0 = NaN){
  if(is.nan(p0)[1,]) {
    Kvector <- kernel_vec(x, data, N, dim, kernel_sigma)
    Xvector <- matrix(rep(x,N),nrow = dim)-t(data)
    GX <- Xvector*matrix(rep(Kvector,dim), nrow = dim, byrow = T)
    p <- mean(Kvector)  #is p a vector or a number
    
    if (adaptive == T){
      kernel_sigma <- kernel_sigma/(p+0.1)
    }
    Kvector <- kernel_vec(x, data, N, dim, kernel_sigma)
    Xvector <- matrix(rep(x,N),nrow = dim)-t(data)
    GX <- Xvector*matrix(rep(Kvector,dim), nrow = dim, byrow = T)
    p <- mean(Kvector)
    g <- -rowMeans(GX)/kernel_sigma^2
    H <- 1/(N*(kernel_sigma^4)) * (GX%*%t(Xvector))-diag(dim)*p/(kernel_sigma^2)
    SI <- -(1/p)*H+(1/p)^2*g%*%t(g)
    
  }
  else {
    
    #the mothord is to compute p on every data poin x_i,then regard vector p as input, comput a new density by 
    #using the equation: p = 1/n*\sum 1/h_j*k(x-x_j/h_j),for example, h_j=h/sqrt(p_j)
    h <- kernel_sigma*(1/sqrt(p0+0.1))
    Kvector <- kernel_vec(x, data, N, dim, kernel_sigma,h = h)
    Xvector <- matrix(rep(x,N),nrow = dim)-t(data)
    GX <- Xvector*matrix(rep(Kvector,dim), nrow = dim, byrow = T)
    p <- as.numeric(Kvector %*% p0)
    g <- -GX %*% p0^2
    H <- 1/(N*(kernel_sigma^4)) * (GX%*%t(Xvector))-diag(dim)*p/(kernel_sigma^2)
    SI <- -(1/p)*H+(1/p)^2*g%*%t(g)
  }
  return(list(p=p, g=g, H=H, SI=SI))
}

#-----------------kernel_vec---------------------------
# Compute the "G_{Sigma_{i}}" vector on point x.
#------------------------------------------------------
kernel_vec <- function(x, data, N, dim, kernel_sigma, h = NaN){
  
  if (is.nan(h)[1,]){
    C <- 1/((2*pi)^(dim/2)*kernel_sigma^dim)
    dx <- (matrix(rep(x,N),nrow = dim)-t(data))/matrix(rep(kernel_sigma,dim*N),nrow = dim, ncol = N)
    
  }
  else {
    C <- 1/((2*pi)^(dim/2)*det(diag(c(1,2,3))))
    dx <- (matrix(rep(x,N),nrow = dim)-t(data))/matrix(rep(h,dim),nrow = dim, ncol = N)
    
  }
  K <- C*exp(-0.5*colSums(dx^2)) #not sum up all the colomn, but compute the sum of every colomn
  return(K)
}