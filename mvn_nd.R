################################################################

# finds MLE by Gradient Descent for any dimension

# This is an extremely inefficient approcah, both theoretically and computationally

# only for testing purposes, shouldn't be used in practice

#################################################################



# function to get the MLE for any n


get_MLE_GD = function(mean,S, n_iter = 100, tolr = 1e-10, start = NULL, alpha = 0.00005){
  
  p = length(mean)
  if(is.null(start) == T ) start = rep(0, p)
  # gradient descent
  
  Sinv = solve(S)  # will be extremely inefficient for higher dimensions
  mu = matrix(0, nrow = n_iter, ncol = p)
  mu[1,] = start
  steps = 1
  
  
  while(steps < n_iter){
    
    D = sweep(x, 2, mu[steps, ], "-")
    score = colSums(D)
    mu[steps+1,] = mu[steps,] + alpha * drop(Sinv %*% score)
    
    diffs = sqrt(sum((mu[steps+1,] - mu[steps,])^2))
    if(diffs < tolr) break
    
    steps = steps + 1
    
  }
  
  print(paste('number of iterations: ', steps))
  
  print('final estimate: ')
  print(mu[steps,])
  
  
  results = list()
  results$n_iter = steps
  results$final = mu[steps,]
  results$estimates = mu
  
  
  
  return(results)
  
}


#--------------
# Example run 
#-------------- 

dims = 3
mu = rep(0,dims) #suitable mean to be given

# to create the covariance matrix
# very specific
# better implementation needed
make_vcov_mat = function(dims){
  s = runif(dims^2,-1,1)
  S = matrix(s, nrow = dims, ncol = dims)
  S =  t(S)%*%S
  S = S/max(eigen(S)$value + 1)
}

S = make_vcov_mat(dims)

# generate samples
x = MASS::mvrnorm(500,mu,S)

# the stepsize(alpha) shoul be very small for this algo
r = get_MLE(mu, S,n_iter = 1000, alpha = 0.00005)


