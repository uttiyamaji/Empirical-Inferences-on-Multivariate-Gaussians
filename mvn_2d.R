# function to draw the contour using the parametrs given and find the MLE of the mean using Gradient Descent


get_MLE = function(mu = c(0,0), rho = 0, s1= 1, s2= 1, n_obs = 500,
                        n_iter = 100, tolr = 1e-10, start = c(0,0), alpha = 0.005){
  
  # generate samples 
  s12 = rho*sqrt(s1)*sqrt(s2)
  S = rbind(c(s1,s12),c(s12,s2))
  x = MASS::mvrnorm(500,mu,S)

  # calculate likelihood and plot the contour over the range
  nx = 40
  ny = 40
  xg = seq(mu[1]-5,mu[1]+5, len = nx)
  yg = seq(mu[2]-5,mu[2]+5, len = ny)
  
  g = expand.grid(xg,yg)
  
  nloglike = function(mu){
    dmv = mvtnorm::dmvnorm(x, mu, S, log = T)
    -sum(dmv)
  }
  
  nLL = apply(g, 1,nloglike)
  
  z = matrix(nLL, nx, ny)
  
  contour(xg, yg, z, nlevels = 40, xlab = expression(mu[1]), ylab = expression(mu[2]))
  abline(h = mu[2], v = mu[1], lty = 2)
  




  # gradient descent

  Sinv = solve(S)  
  mu = matrix(0, nrow = n_iter, ncol = 2)
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
  
 points(mu, pch = 20)
 
 return(results)
  
}









