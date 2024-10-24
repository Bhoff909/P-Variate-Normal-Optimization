library(matrixcalc)

sqrtm <- function (A) {
  # Obtain matrix square root of a matrix A
  a = eigen(A)
  sqm = a$vectors %*% diag(sqrt(a$values)) %*% t(a$vectors)
  sqm = (sqm+t(sqm))/2
}
# a
gen <- function(n,p,mu,sig,seed = 2023){
  #---- Generate data from a p-variate normal with mean mu and covariance sigma
  # mu should be a p by 1 vector
  # sigma should be a positive definite p by p matrix
  # Seed can be optionally set for the random number generator
  set.seed(seed)
  # generate data from normal mu sigma
  z = matrix(rnorm(n*p),n,p)
  datan = z %*% sqrtm(sig) + matrix(mu,n,p, byrow = TRUE)
  datan
}

#Write separate functions to compute the log-likelihood, the gradient of the log-likelihood, 
#and the Hessian of the log-likelihood. Also write separate functions for each of the methods: 
#steepest ascent, Newton, and Fisher-scoring

#should sums be computed outside of functions?

log_likelihood <- function(data,mu,sig,n,p){
  siginv = solve(sig)
  C= matrix(0,p,p); # initializing sum of (xi-mu)(xi-mu)^T
  for (i in 1:n){
    xm = data[i,] - mu
    C = C + xm %*% t(xm)
  }
  if(is.null(siginv)){
    log_det_sig <- log(det(sig))
  } else {
    #-- in this case siginv is input so we use the fact that det(sig)=1/det(siginv)
    log_det_sig <-log(1/det(siginv))
  }
  l = -(n*p*log(2*pi)+n*log_det_sig + sum(siginv * C ))/2
}

gradient <- function(data, mu, sig){
  p = dim(data)[2]
  siginv = solve(sig)
  
  C= matrix(0,p,p); # initializing sum of (xi-mu)(xi-mu)^T
  sxm = matrix(0,p,1) # initializing sum of xi-mu
  gradmu = sxm; # initializing this sum is used for the gradient w.r.t. mu
  for (i in 1:n){
    xm = data[i,] - mu
    sxm = sxm + xm
    C = C + xm %*% t(xm)
  }
  gradmu = siginv %*% sxm
  gradsig = n*siginv - siginv %*% C %*% siginv
  
    for (i in 1:p){
      for (j in 1:p){
        if(i == j){gradsig[i,j] = -1/2*gradsig[i,j]}else{gradsig[i,j] = -1*gradsig[i,j]}
      }
    }
  #ensure symmetry of sigma
  gradsig = vec2mat(mat2vec(gradmu, gradsig, p),p)$sigma
  
  gradlist <- list("gradmu" = gradmu, "gradsig" = gradsig)
  return(gradlist)
}

n <- 200
p <- 3
sig <- matrix(c(1,.7,.7,.7,1,.7,.7,.7,1),3,3) # known sigma Note <<- makes it global
mu <- matrix(c(-1,1,2))
datan = gen(n,p,mu,sig)

mu0 <- matrix(c(0,0,0))
sig0 <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)

testout = gradient(datan, mu0, sig0)
testout$gradmu
testout$gradsig

!is.positive.definite(testout$gradsig)

mu1 = mu + testout$gradmu
sig1 = sig + testout$gradsig

mu1
sig1
halve = 0
while(!is.positive.definite(sig1)){
  halve = halve + 1
  mu1 = mu + (testout$gradmu)/2^halve  # Steepest Ascent
  sig1 = sig + (testout$gradsig)/2^halve}
mu1
sig1

log_likelihood(datan, mu1, sig1,n,p)

gradvec = mat2vec(testout$gradmu, testout$gradsig, 3)
gradvec

hessian <- function(data, mu, sig){
  p = dim(data)[2]
  n = dim(data)[1]

  siginv = solve(sig)
  I = diag(p)
  
  s= matrix(0,p,p); # initializing sum of (xi-mu)(xi-mu)^T
  sxm = matrix(0,p,1) # initializing sum of xi-mu
  
  for (i in 1:n){
    xm = data[i,] - mu
    sxm = sxm + xm
    s = s + (xm %*% t(xm))
  }
  c = siginv %*% sxm
  z = (-n*I + 2*siginv %*% s) %*% siginv
  print(z)
  stop()

  mumu = -n*siginv
  
  columnindex = 1
  rowindex = 1
  sigsig = matrix(0, p*(p+1)/2, p*(p+1)/2)
  for(i in 1:p){
    for(j in 1:i){
      for(k in 1:p){
        for(l in 1:k){
          if(i == j & l == k){sigsig[rowindex, columnindex] = -1/2*(z[k,i]*siginv[i,k])}
          if(i != j & k != l){sigsig[rowindex, columnindex] = -1/2*(z[k,i]*siginv[j,l] + z[l,j]*siginv[i,k] + z[k,j]*siginv[i,l] + z[l,i]*siginv[j,k])}
          if(i != j & k == l){sigsig[rowindex, columnindex] = -1/2*(z[k,i]*siginv[j,k] + z[k,j]*siginv[i,k])}
          if(i == j & k != l){sigsig[rowindex, columnindex] = -1/2*(z[l,i]*siginv[i,k] + z[k,i]*siginv[i,l])}
          if(rowindex < p*(p+1)/2){rowindex = rowindex + 1}else{rowindex = 1}
          }
      }
      if(columnindex < p*(p+1)/2){columnindex = columnindex + 1}else{columnindex = 1}
    }
  }

  
  sigmu = matrix(0, p, p*(p+1)/2)
  columnindex = 1
  #new matrix for mu sigma section
  for(l in 1:p){
    for(k in 1:l){
      for(i in 1:p){
        if(l == k){sigmu[i,columnindex] = -siginv[i,k]*c[k]}else{sigmu[i,columnindex] = -siginv[i,k]*c[l] - siginv[i,l]*c[k]}
      }
      if(columnindex < p*(p+1)/2){columnindex = columnindex + 1}else{columnindex = 1}
    }
  }
  

  sigmuleft = t(sigmu)
  
  combination = cbind(mumu, sigmu)
  combination2 = cbind(sigmuleft, sigsig)
  finalcombo = rbind(combination, combination2)
  
  return(finalcombo)
}

######################test hessian
print(hessian(datan, mutest, sigtest))

sigtest <- matrix(c(1,.5,.5,.5,1,.5,.5,.5,1),3,3) # known sigma Note <<- makes it global
mutest <- matrix(c(-1.5,1.5,2.3))
######################test hessian
print(hessian(datan, mutest, sigtest))

mat2vec <- function(mu, sig, p){
  #vectorize mu
  out = c(mu)
  
  #vectorize lower diagonals of sigma
  sigvec = c(sig[1:1])
  for (i in 2:p){
    sigvec = append(sigvec, sig[i, 1:i])
  }
  
  out = append(out, sigvec)
  return(out)
}

vec = mat2vec(mu, sig, 3)
vec
mu2 = matrix(vec[1:p])

#######test calculating second derivatives************
#for sigma
sigmu = matrix(0, p, p*(p+1)/2)
columnindex = 1
#new matrix for mu sigma section
for(i in 1:p){
  for(k in 1:i){
    for(l in 1:p){
      #if(k == l){sigmu[k,columnindex] = siginv[i,k]*c[k]}else{sigmu[k,columnindex] = siginv[i,k]*c[l] + siginv[i,l]*c[k]}
      if(i == k){sigmu[l,columnindex] = siginv[i,k]*c[i]}else{sigmu[l,columnindex] = siginv[i,k]*c[k] + siginv[i,l]*c[i]}
      #if(i == k){sigmu[l,columnindex] = 1}else{sigmu[l,columnindex] = 2}
    }
    if(columnindex < 6){columnindex = columnindex + 1}else{columnindex = 1}
  }
}
sigmu
columnindex = 1
for(k in 1:p){
  for(l in 1:i){
    for(i in 1:p){
      #if(k == l){sigmu[k,columnindex] = siginv[i,k]*c[k]}else{sigmu[k,columnindex] = siginv[i,k]*c[l] + siginv[i,l]*c[k]}
      if(k==l){sigmu[i,columnindex] = 1}else{sigmu[i,columnindex] = 2}
    }
    if(columnindex < 6){columnindex = columnindex + 1}else{columnindex = 1}
  }
}
sigmu

lowerleft = matrix(0, p*(p+1)/2,p)
rowindex = 1
for(i in 1:p){
  for(j in 1:i){
    for(k in 1:p){
      if(i == j){lowerleft[rowindex,k] = 1}else{lowerleft[rowindex,k] = 2}
    }
    if(rowindex < 6){rowindex = rowindex + 1}else{rowindex = 1}
  }
}
lowerleft

sigtest <- matrix(c(1,.5,.5,.5,1,.5,.5,.5,1),3,3) # known sigma Note <<- makes it global

mutest <- matrix(c(-1.5,1.5,2.3))
#####################################################


vec2mat <- function(thetvec, p){
  mumat = matrix(thetvec[1:p])
  sigmat = matrix(0, nrow = p, ncol = p)
  
  pos = p+1
  for (i in 1:p){
    sigmat[i, 1:i] = thetvec[c(pos:(pos+i-1))]
    pos = pos + i
  }
  #return upper diagonals
  for(i in 1:p){
    for(j in i : p){
      if(j > i){sigmat[i,j] = sigmat[j,i]}
    }
  }

  matlist <- list("mu" = mumat, "sigma" = sigmat)
  return(matlist)
}

outmat = vec2mat(vec, p)
sig2 = outmat$sigma
sig2
mu2 = outmat$mu
#--- generating data and graphing
# Generate data inside f, (mu1 and mu2 can be entered as vectors, 
# if evaluation at several values of mu1 and mu2 is to be done)
n <- 200
p <- 3
sig <<- matrix(c(1,.7,.7,.7,1,.7,.7,.7,1),3,3) # known sigma Note <<- makes it global
siginv <- solve(sig)
mu <- matrix(c(-1,1,2))
datan = gen(n,p,mu,sig)

mu0 <- matrix(c(0,0,0))
sig0 <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)


SA <- function (mu, datan, sig, maxit, tolerr, tolgrad) {
  n <- nrow(datan)
  p = dim(datan)[2]
  
  for (it in 1:maxit){
    a <- log_likelihood(datan, mu, sig,n,p)
    agrad <- gradient(datan,mu,sig)
    
    mu1 = mu + agrad$gradmu
    sig1 = sig + agrad$gradsig
    halve = 0;
    
    print("iteration  halving    log-likelihood   ||gradient||")
    
    #check if sig1 is positive definite before sending to log likelihood
    while(!is.positive.definite(sig1) & halve <= 20){
      if(it == 1 || it == 2 || it == 425 || it == 426){print(sprintf('%2.0f          %2.0f         ',it, halve))}
      halve = halve + 1
      mu1 = mu + (agrad$gradmu)/2^halve  # Steepest Ascent
      sig1 = sig + (agrad$gradsig)/2^halve
    }
    if (halve >= 20) print('Step-halving failed after 20 halvings')
    
    atmp = log_likelihood(datan, mu1, sig1,n,p)
    
    #continue halving if new value smaller than prior
    while (atmp < a & halve <= 20){
      #L2 norm of gradient
      gradnorm <- norm(mat2vec(gradient(datan, mu1, sig1)$gradmu,gradient(datan, mu1, sig1)$gradsig,p),type = "2")
      print(sprintf('%2.0f          %2.0f         %2.4f         %0.1e',it, halve,atmp,gradnorm))
      halve = halve+1
      mu1 = mu + (agrad$gradmu)/2^halve  # Steepest Ascent
      sig1 = sig + (agrad$gradsig)/2^halve
      atmp = log_likelihood(datan, mu1, sig1,n,p)
    }
    if (halve >= 20) print('Step-halving failed after 20 halvings')

    #modified relative error
    mre = max(((mat2vec(mu1, sig1, p)) - (mat2vec(mu, sig, p)))/pmax(1,abs(mat2vec(mu1, sig1, p))))
    #L2 norm of gradient
    gradnorm = norm(mat2vec(gradient(datan, mu1, sig1)$gradmu,gradient(datan, mu1, sig1)$gradsig,p),type = "2")
    
    print(sprintf('%2.0f          %2.0f         %2.4f         %0.1e',it, halve,atmp,gradnorm))
    print('-----------------------------')
    if(mre < tolerr & gradnorm < tolgrad){break}
    
    mu = mu1
    sig = sig1

  }
  paramlist <- list("mu" = mu, "sigma" = sig)
  return(paramlist)
}

#b
#initial values
mu0 <- matrix(c(0,0,0))
sig0 <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)
SA(mu0, datan, sig0, 500, 1e-6, 1e-5)

newton <-function(mu, datan, sig, maxit, tolerr, tolgrad){
    
   n = dim(datan)[1]
   p = dim(datan)[2]
  
   for(it in 1:maxit){
    a <- log_likelihood(datan, mu, sig,n,p)
    grad = gradient(datan, mu, sig)
    gradvec = mat2vec(grad$gradmu, grad$gradsig, p)
    dir = -1*solve(hessian(datan, mu, sig)) %*% gradvec
    dirmat = vec2mat(dir, p)
    
    mu1 = mu + dirmat$mu
    sig1 = sig + dirmat$sigma
    
    halve = 0;
    
    print("iteration  halving    log-likelihood   ||gradient||")
    #check if sig1 is positive definite before sending to log likelihood
    while(!is.positive.definite(sig1) & halve <= 20){
      halve = halve + 1
      mu1 = mu + (dirmat$mu)/2^halve
      sig1 = sig + (dirmat$sigma)/2^halve
      print(sprintf('%2.0f          %2.0f           %2.4f', a, it, halve))
    }
    
    atmp = log_likelihood(datan, mu1, sig1,n,p)
    
    #continue halving if new value smaller than prior
    while (atmp < a & halve <= 20){
      halve = halve+1
      mu1 = mu + (dirmat$mu)/2^halve
      sig1 = sig + (dirmat$sigma)/2^halve
      atmp = log_likelihood(datan, mu1, sig1,n,p)
      #L2 norm of gradient
      gradnorm <- norm(mat2vec(gradient(datan, mu1, sig1)$gradmu,gradient(datan, mu1, sig1)$gradsig,p),type = "2")
      print(sprintf('%2.0f          %2.0f         %2.4f         %0.1e',it, halve,atmp,gradnorm))
    }
    if (halve >= 20) print('Step-halving failed after 20 halvings')
   
     #modified relative error
    mre = max(((mat2vec(mu1, sig1, p)) - (mat2vec(mu, sig, p)))/pmax(1,abs(mat2vec(mu1, sig1, p))))
    #L2 norm of gradient
    gradnorm = norm(mat2vec(gradient(datan, mu1, sig1)$gradmu,gradient(datan, mu1, sig1)$gradsig,p),type = "2")
    if(mre < tolerr & gradnorm < tolgrad){break}
    
    print(sprintf('%2.0f          %2.0f         %2.4f         %0.1e',it, halve,atmp,gradnorm))
    print('-----------------------------')
    
    mu = mu1
    sig = sig1

   }
   paramlist <- list("mu" = mu, "sigma" = sig)
   return(paramlist)
}

newton(mu0, datan, sig0, 500, 1e-7, 1e-7)

fisher_info <-function(data, sig){
  n = dim(data)[1]
  print(n)
  siginv = solve(sig)
  topleft = -n*siginv
  
  topright = matrix(0,3,6)
  botleft = t(topright)

  columnindex = 1
  rowindex = 1
  botright = matrix(0, p*(p+1)/2, p*(p+1)/2)
  for(i in 1:p){
    for(j in 1:i){
      for(k in 1:p){
        for(l in 1:k){
          if(i == j & l == k){botright[rowindex, columnindex] = -n/2*(siginv[k,i]*siginv[i,k])}
          if(i != j & k != l){botright[rowindex, columnindex] = -n/2*(siginv[k,i]*siginv[j,l] + siginv[l,j]*siginv[i,k] + siginv[k,j]*siginv[i,l] + siginv[l,i]*siginv[j,k])}
          if(i != j & k == l){botright[rowindex, columnindex] = -n/2*(siginv[k,i]*siginv[j,k] + siginv[k,j]*siginv[i,k])}
          if(i == j & k != l){botright[rowindex, columnindex] = -n/2*(siginv[l,i]*siginv[i,k] + siginv[k,i]*siginv[i,l])}
          if(rowindex < p*(p+1)/2){rowindex = rowindex + 1}else{rowindex = 1}
        }
      }
      if(columnindex < p*(p+1)/2){columnindex = columnindex + 1}else{columnindex = 1}
    }
  }
  
  top = cbind(topleft, topright)
  bottom = cbind(botleft, botright)
  finalcombo = rbind(top, bottom)
  return(finalcombo)
  
}
sigtest
solve(fisher_info(datan, sigtest))

fisher <-function(mu, datan, sig, maxit, tolerr, tolgrad){
  n <- nrow(datan)
  p = dim(datan)[2]
  
  for(it in 1:maxit){
    a <- log_likelihood(datan, mu, sig,n,p)
    grad = gradient(datan, mu, sig)
    gradvec = mat2vec(grad$gradmu, grad$gradsig, p)
    info = fisher_info(datan, sig)
    dir = -solve(info) %*% gradvec
    dirmat = vec2mat(dir, p)
    
    mu1 = mu + dirmat$mu
    sig1 = sig + dirmat$sigma
    
    halve = 0;
    
    #check if sig1 is positive definite before sending to log likelihood
    while(!is.positive.definite(sig1) & halve <= 20){
      print(sprintf('%2.0f          %2.0f         %2.4f',it, halve,a))
      halve = halve + 1
      mu1 = mu + (dirmat$mu)/2^halve
      sig1 = sig + (dirmat$sigma)/2^halve
    }
    
    atmp = log_likelihood(datan, mu1, sig1,n,p)
    
    #continue halving if new value smaller than prior
    while (atmp < a & halve <= 20){
      #L2 norm of gradient
      gradnorm <- norm(mat2vec(gradient(datan, mu1, sig1)$gradmu,gradient(datan, mu1, sig1)$gradsig,p),type = "2")
      print(sprintf('%2.0f          %2.0f         %2.4f         %0.1e',it, halve,atmp,gradnorm))
      halve = halve+1
      mu1 = mu + (dirmat$mu)/2^halve
      sig1 = sig + (dirmat$sigma)/2^halve
      atmp = log_likelihood(datan, mu1, sig1,n,p)
    }
    if (halve >= 20) print('Step-halving failed after 20 halvings')
    
    #modified relative error
    mre = max(((mat2vec(mu1, sig1, p)) - (mat2vec(mu, sig, p)))/pmax(1,abs(mat2vec(mu1, sig1, p))))
    #L2 norm of gradient
    gradnorm = norm(mat2vec(gradient(datan, mu1, sig1)$gradmu,gradient(datan, mu1, sig1)$gradsig,p),type = "2")
    if(mre < tolerr & gradnorm < tolgrad){break}
    
    print("iteration  halving    log-likelihood   ||gradient||")
    print(sprintf('%2.0f          %2.0f         %2.4f         %0.1e',it, halve,atmp,gradnorm))
    print('-----------------------------')
    
    
    mu = mu1
    sig = sig1

  }
  paramlist <- list("mu" = mu, "sigma" = sig)
  return(paramlist)
}

fisher(mu0, datan, sig0, 500, 1e-7, 1e-7)

sig0 <- matrix(c(1,.5,.5,.5,1,.5,.5,.5,1),3,3) # known sigma Note <<- makes it global
mu0 <- matrix(c(-1.5,1.5,2.3))

