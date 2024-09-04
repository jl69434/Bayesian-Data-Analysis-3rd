# Replicating Figure 13.7 in BDA3 (Expectation propagation)


num = 4
iter = 10
k = 1

x = matrix(c(1, -0.86, 1, -0.30, 1, -0.05, 1, 0.73), 2, 4)
m = c(5, 5, 5, 5)
y = c(0, 1, 3, 5)

cinvmu.app = cinvmu.cav = array(0, c(2, 1, iter+1, num))
cinv.app = cinv.cav = array(0, c(2, 2, iter+1, num))
cinvmu = array(0, c(2, 1, num*iter+1))
cinv = array(0, c(2, 2, num*iter+1))

# Set initial value to get the same result in BDA3
cinv[ , ,1] = c(4, 0, 0, 4)
cinv.app[ , , 1, ] = c(1, 0, 0, 1)


for(j in 1: iter){
  for(i in 1:num){
    # Step 1
    cinvmu.cav[ , , j+1, i] = cinvmu[ , , k] - cinvmu.app[ , , j, i]
    cinv.cav[ , ,j+1, i] =  cinv[ , , k] - cinv.app[ , , j, i]
    
    # Step 2
    M.cav = as.numeric(t(x[,i]) %*% solve(cinv.cav[ , , j+1, i]) %*% cinvmu.cav[ , , j+1, i])
    V.cav = as.numeric(t(x[,i]) %*% solve(cinv.cav[ , , j+1, i]) %*% x[,i])
    
    # Step 3
    tilted.0 = function(eta){
      dnorm(eta, M.cav, sqrt(V.cav)) * 
        factorial(m[i]) / (factorial(y[i]) * factorial(m[i]-y[i])) * 
        (exp(eta)/(1+exp(eta)))^y[i] * (1/(1+exp(eta)))^(m[i]-y[i])
    }
    tilted.1 = function(eta){
      eta * dnorm(eta, M.cav, sqrt(V.cav)) * 
        factorial(m[i]) / (factorial(y[i]) * factorial(m[i]-y[i])) * 
        (exp(eta)/(1+exp(eta)))^y[i] * (1/(1+exp(eta)))^(m[i]-y[i])
    }
    tilted.2 = function(eta){
      eta^2 * dnorm(eta, M.cav, sqrt(V.cav)) * 
        factorial(m[i]) / (factorial(y[i]) * factorial(m[i]-y[i])) * 
        (exp(eta)/(1+exp(eta)))^y[i] * (1/(1+exp(eta)))^(m[i]-y[i])
    }
    
    E0 =  integrate( tilted.0, M.cav-10*sqrt(V.cav), M.cav+10*sqrt(V.cav) )$value
    E1 =  integrate( tilted.1, M.cav-10*sqrt(V.cav), M.cav+10*sqrt(V.cav) )$value
    E2 =  integrate( tilted.2, M.cav-10*sqrt(V.cav), M.cav+10*sqrt(V.cav) )$value
    M = E1/E0 ; V = E2/E0 - (E1/E0)^2
    
    # Step 4
    V.app = 1 / (1/V - 1/V.cav)
    M.app = V.app * (M/V - M.cav/V.cav)
    
    # Step 5 (update for mu and covariance for i)
    cinv.app[ , , j+1, i] = (1/V.app) * x[,i] %*% t(x[,i])
    cinvmu.app[ , , j+1, i] = (M.app/V.app) * x[,i] 
    
    # Step 6 (update mu and covariance)
    k = k+1
    cinvmu[ , , k] = cinvmu.cav[ , , j+1, i] + cinvmu.app[ , , j+1, i]
    cinv[ , , k] = cinv.cav[ , , j+1, i] + cinv.app[ , , j+1, i]
  }
}

sigma1 = sigma2 = mu1 = mu2 = corr = matrix(0, 41, 1)
cov.matrix = array(0, c(2, 2, 41))
for(i in 1:41){
  cov.matrix[ , , i] = solve(cinv[ , , i])
  sigma1[i] = sqrt(cov.matrix[1,1,i])
  sigma2[i] = sqrt(cov.matrix[2,2,i])
  mu1[i] = (cov.matrix[ , ,i] %*% cinvmu[ , , i])[1,1]
  mu2[i] = (cov.matrix[ , ,i] %*% cinvmu[ , , i])[2,1]
  corr[i] = cov.matrix[1,2,i]/( sqrt(cov.matrix[1,1,i]) * sqrt(cov.matrix[2,2,i]))
}


# Plotting Figure 13.7
dev.off()
par(mfrow=c(2,3))
plot(mu1, type='l', xlab='', xaxt='n', ylab='',  yaxt='n', ylim=c(-1.1, 1.8), bty='l', main='Updating of ¥ì1')
axis(side=1, at=c(0, 20, 40), labels=c(0, 5, 10))
axis(side=2, at=c(-1, 0, 1))
plot(sigma1, type='l', xlab='', ylab='', xaxt='n', yaxt='n', ylim=c(0, 1.2), bty='l', main='Updating of ¥ò1')
axis(side=1, at=c(0, 20, 40), labels=c(0, 5, 10))
axis(side=2, at=c(0, 0.5, 1))
plot.new()
plot(mu2, type='l', xlab='iteration', xaxt='n', ylab='', yaxt='n', ylim=c(-1, 12), bty='l', main='Updating of ¥ì2')
axis(side=1, at=c(0, 20, 40), labels=c(0, 5, 10))
axis(side=2, at=c(0, 5, 10))
plot(sigma2, type='l', xlab='iteration', xaxt='n', ylab='', yaxt='n', ylim=c(0,6), bty='l', main='Updating of ¥ò2')
axis(side=1, at=c(0, 20, 40), labels=c(0, 5, 10))
axis(side=2, at=c(0, 2.5, 5))
plot(corr, type='l', xlab='iteration', xaxt='n', ylab='', yaxt='n', ylim=c(0,0.7), bty='l', main='Updating of ¥ñ')
axis(side=1, at=c(0, 20, 40), labels=c(0, 5, 10))
axis(side=2, at=c(0, 0.25, 0.5))