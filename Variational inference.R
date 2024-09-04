# Replicating Figure 13.5 and 13.6 in BDA3 (the progress of the variational algorithm)

y = c(28, 8, -3, 7, -1, 1, 18, 12)
sigma = c(15, 10, 16, 11, 9, 11, 10, 18)
sig.sq = sigma^2

iter = 100
M = S = matrix(0, iter, 9)
M.tau = matrix(0, iter, 1)
KL = matrix(0, iter, 1)

# Set initial value to get the same result in BDA3
M[1, ] = rnorm(1)
S[1, ] = runif(1)^2

# According to BDA3, we can set initial value as follows
#for(i in 1:9){
#  M[1,i] = rnorm(1)
#  S[1,i] = runif(1)^2
#}

# Calculate initial values
M.tau[1] = 1/7*sum((M[1,1:8] - M[1,9])^2) + 1/7*sum(S[1,1:8] + S[1,9])
KL[1] = 0.5*sum(((y[1:8] - M[1,1:8])^2 + S[1,1:8]) / sig.sq[1:8])
       +0.5*sum(((M[1,1:8] - M[1,9])^2 + S[1,1:8] + S[1,9]) / M.tau[1])
       -sum(log(sqrt(S[1,1:9])))


# Updating procedure
for(i in 2:iter){
  for(j in 1:8){
    
  # Step (13.18)
    M[i,j] = (1/sig.sq[j] * y[j] + 1/M.tau[i-1]*M[i-1,9]) / (1/sig.sq[j] + 1/M.tau[i-1])
    S[i,j] = 1 / (1/sig.sq[j] + 1/M.tau[i-1])
  }
  
  # Step (13.19)
  M[i,9] = 1/8 * sum(M[i, 1:8])
  S[i,9] = 1/8 * M.tau[i-1]
  
  # Step (13.20)
  M.tau[i] = 1/7*sum((M[i,1:8] - M[i,9])^2) + 1/7*sum(S[i,1:8] + S[i,9])
  
  # Calculate KL divergence
  KL[i] = 0.5 * sum(((y[1:8] - M[i,1:8])^2 + S[i,1:8]) / sig.sq[1:8])
         +0.5 * sum(((M[i,1:8] - M[i,9])^2 + S[i,1:8] + S[i,9]) / M.tau[i])
         -sum(log(sqrt(S[i,1:9])))
}



# Plotting Figure 13.5
dev.off()
par(mfrow=c(2,3))
name = LETTERS[1:8]
plot(M[,9], type='l', xlab='', xaxt='n', ylab='', ylim = c(0,10), main='Updating of Mレ', bty='l')
axis(side=1, at=c(0, 50, 100), labels=c("", "", ""))
plot(sqrt(S[,9]), type='l', xlab='', xaxt='n', ylab='', ylim = c(0,10), main='Updating of Sレ', bty='l')
axis(side=1, at=c(0, 50, 100), labels=c("", "", ""))
plot(sqrt(M.tau), type='l', xlab='', xaxt='n', ylab='', ylim = c(0,10), main='Updating of Mン', bty='l')
axis(side=1, at=c(0, 50, 100), labels=c("", "", ""))
plot((M[,1]), type='l', ylim=c(-5,15), xlab='iteration', xaxt='n', ylab='', main='Updating of Mメ', bty='l')
axis(side=1, at=c(0, 50, 100))
for(i in 1:8){
  lines((M[,i]), type='l', ylim=c(-5,15))
  text(x=102, y=M[iter,i], labels=name[i], cex=1)
}
plot(sqrt(S[,1]), type='l', ylim=c(0,10), xlab='iteration', xaxt='n', ylab='', main='Updating of Sメ', bty='l')
axis(side=1, at=c(0, 50, 100))
for(i in 1:8){
  dot = c(101, 101, 101, 101, 101, 103, 103, 101)
  lines(sqrt(S[,i]), type='l', ylim=c(0,10))
  text(x=dot[i], y=sqrt(S[iter,i]), labels=name[i], cex=1)
}
plot(KL, type='l', xlab='iteration', xaxt='n', ylab='', main = '(Unscaled) Kullback-Leibler divergence', bty='l')
axis(side=1, at=c(0, 50, 100))



# Plotting Figure 13.6
dev.off()
par(mfrow=c(1,3))
quantile = c(0.5, 0.25, 0.75, 0.05, 0.95)
x = seq(0, 100, length=iter)
y = array(0, c(iter, 5, 3))
for(i in 1:3){
  for(j in 1:5){
    y[, j, i] = qnorm(quantile[j], mean=M[,i], sd=sqrt(S[,i]))
  }
title = paste("Variational inference for メ", i, sep="")
plot(y[, 1, i], xlab='iteration', xaxt='n', ylab='', type='l', ylim=c(-10,30), bty='l', main = title)
axis(side=1, at=c(0, 50, 100))
polygon(c(x, rev(x)), c(y[, 4, i], rev(y[, 5, i])),  col = gray(0.90), border=NA)
polygon(c(x, rev(x)), c(y[, 2, i], rev(y[, 3, i])),  col = gray(0.65), border=NA)
lines(y[, 1, i], type='l')
}

