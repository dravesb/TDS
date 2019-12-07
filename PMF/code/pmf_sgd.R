#--------------------------------------
#
#       Image Processing
#
#--------------------------------------
library("png")
library("colorspace")

#read in image
img <- readPNG("~/Documents/Work/TDS/PMF/data/pulisic.png")

#check dimension of data frame
dim(img)

#get color dictionary
y <- rgb(img[,,1], img[,,2], img[,,3], alpha = img[,,4])

#make greyscale
yg <- desaturate(y)
yn <- col2rgb(yg)[1, ]

#set image dimensions
dim(y) <- dim(yg) <- dim(yn) <- dim(img)[1:2]

#reshape yn 
yn <- yn[dim(yn)[1]:1,1:dim(yn)[2]]

heatmap(yn,Rowv=NA,Colv=NA,col=paste("gray",1:99,sep=""),
        labRow = FALSE, labCol = FALSE)

#--------------------------------------
#
#       Corrupt Image
#
#--------------------------------------
set.seed(1985)

#sample x's to corrup
alpha <- .3
ind.x <- sample(1:dim(yn)[1], size = round(alpha*dim(yn)[1]))
ind.y <- sample(1:dim(yn)[2], size = round(alpha*dim(yn)[2]))

#corrupt image
yn_corrupt <- yn
yn_corrupt[ind.x, ind.y] <- NA

#plot image
heatmap(yn_corrupt,Rowv=NA,Colv=NA,col=paste("gray",1:99,sep=""),
        labRow = FALSE, labCol = FALSE)

#--------------------------------------
#
#           PMF Code 
#
#--------------------------------------

#set up gradient functions
grad_u <- function(X, V, u, i, sigma2, sigma2_u){
  
  #get the full sum (may want to sample for stochastic gradients)
  sum1 <- numeric(d)
  
  for(j in 1:(nrow(X))){
    if(I[i,j] == 1){ #is the value missing? 
      sum1 <- sum1 + as.numeric(X[i,j] - crossprod(u, V[j,])) * V[j,]  
    }
  }
  
  #return the gradient
  return(-sum1/sigma2 + u/sigma2_u)
  
}
grad_v <- function(X, v, U, j, sigma2, sigma2_v){
  
  #get the full sum (may want to sample for stochastic gradients)
  sum1 <- 0 
  for(i in 1:nrow(U)){
    if(I[i,j] == 0){ #is the value missing? 
      sum1 <- sum1 + 0
    }else{
      sum1 <- sum1 + as.numeric(X[i,j] - crossprod(U[i,], v)) * U[i,]  
    }
  }
  
  #return the gradient
  return(-sum1/sigma2 + v/sigma2_v)
  
}

#probability matrix factorization function 
pmf <- function(X, U, V, eta = 1e-5){
  
  #initialize variables
  d <- ncol(U)
  sigma2 <- 1
  sigma2_u <- 1
  sigma2_v <- 1
  I <- matrix(1, ncol = ncol(X), nrow = nrow(X))
  I[which(is.na(X), arr.ind = TRUE)] <- 0
  max.iters <- 10000
  Xhat <- tcrossprod(U, V)
  
  #convergence criterion
  iter <- 1
  conv.crit <- FALSE
  
  while(!conv.crit){
    
    #update U rows
    for(i in 1:nrow(U)){
      U[i,] <- U[i,] - eta * grad_u(X, V, U[i,], i, sigma2, sigma2_u)
    }
    
    #update V rows
    for(j in 1:nrow(V)){
      V[j,] <- V[j,] - eta * grad_v(X, V[j,], U, j, sigma2, sigma2_v)
    }
    
    #update convergence criterion
    Xhat.new <- tcrossprod(U, V)
    
    par(mfrow = c(3, 1))
    hist(Xhat.new, main = paste("Iter:", iter)); hist(U, main = "U"); hist(V, main = "V")
    
    conv.crit <- ifelse(iter > max.iters || norm(Xhat - Xhat.new, type = "F") < 1e-5, TRUE, FALSE)
    
    print(norm(Xhat - Xhat.new, type = "F"))
    
    #update parameters
    Xhat <- Xhat.new
    iter <- iter +1 
    print()
    
    
  }
  
  #return the two estimated matrices
  return(list(U, V))
  
}

#--------------------------------------
#
#           Get model estimates 
#
#--------------------------------------

#determine embedding dimension d 
plot(svd(yn)$d)

#let's try d = 5 
d <- 5

#set initial values to average greyscale values
U <- matrix(sqrt(mean(yn)/d), nrow = nrow(yn), ncol = d)
V <- matrix(sqrt(mean(yn)/d), nrow = ncol(yn), ncol = d)

res <- pmf(yn_corrupt, U, V, eta = 1e-5)





