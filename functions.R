require(matrixStats)

readXyw <- function() {
  #Reads X,y,w of given clad problem
  require(pracma)
  #X <- read.table("X.txt")
  #X <- read.table("X_200_4.txt")
  #X <- read.table("X_200_5.txt")
  X <- read.table("X_601_5.txt")
  X <- X[,-1]
  X <- as.matrix(X)
  #y <- read.table("ys.txt")
  #y <- read.table("ys_200_4.txt")
  #y <- read.table("ys_200_5.txt")
  y <- read.table("ys_601_5.txt")
  y <- y[,-1]
  y <- as.vector(y)
  w <- repmat(1,dim(X)[1],1)
  results1 <- list(X = X, y = y, w = w)
  return(results1)  
}


standardizeX <- function(X) {
  require(pracma)
  #Standardizes X
  mu <- colMeans(X)
  sigma <- colSds(X)
  testX <- (X - repmat(mu,dim(X)[1],1)) / repmat(sigma,dim(X)[1],1) 
  p <- dim(X)[2]
  for (j in 1:p) {
    if (is.nan(testX[1,j])) {
      X[,j] <- X[,j]
    }
    else {
      X[,j] <- testX[,j]
    }
  }
  results2 <- list(X = X, mu = mu, sigma = sigma)
  return(results2)  
}


definecAb <- function(X,y,w) {
  
  #Defines c,A,b for milp.R
  require(pracma)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  c3 <- repmat(0,1,n)  # gammas integer
  c1 <- repmat(0,1,p)    #betas
  c2 <- repmat(0,1,n)      #phis
  c4 <- repmat(1,1,n)        #sm's
  c5 <- repmat(1,1,n)          #sp's
  
  c<- cbind(c3,c1,c2,c4,c5)
  

  d=50
  #d=40
  #d=20
  #d=10
  #d=5
  M <- numeric(n)
  for (i in 1:n){
    M[i] <- abs(X[i,1])*d+abs(X[i,-1])%*%t(repmat(d,1,p-1))
  }

  A1 <- matrix(0,n,4*n+p)
  for (i in 1:n) {    #constraint 3c in pdf ssrn
    A1[i,1:n] <- 0   #gammas
    for (j in 1:p) {
      A1[i,(n+j)] <- X[i,j]   #betas
    }
    A1[i,(n+p+1):(2*n+p)] <- -1*(1:n==i)   #phis
    A1[i,(2*n+p+1):(3*n+p)] <- 0          #sm's
    A1[i,(3*n+p+1):(4*n+p)] <- 0           #sp's
  }
  
  b1 <- numeric(n)
  
  A2 <- matrix(0,n,4*n+p)
  for (i in 1:n) {    #constraint 3e in pdf ssrn
    A2[i,1:n] <- M[i]*(1:n==i)   #gammas
    for (j in 1:p) {
      A2[i,(n+j)] <- -X[i,j]   #betas
    }
    A2[i,(n+p+1):(2*n+p)] <- 1*(1:n==i)   #phis
    A2[i,(2*n+p+1):(3*n+p)] <- 0   #sm's
    A2[i,(3*n+p+1):(4*n+p)] <- 0    #sp's
  }
  
  b2 <- numeric(n)
  for (i in 1:n) {
    b2[i] <- M[i]
  }
  
  A3 <- matrix(0,n,4*n+p)
  for (i in 1:n) {   #constraint 3f in pdf ssrn
    A3[i,1:n] <- -M[i]*(1:n==i)   #gammas
    A3[i,(n+1):(n+p)] <- 0   #betas
    A3[i,(n+p+1):(2*n+p)] <- 1*(1:n==i)   #phis
    A3[i,(2*n+p+1):(3*n+p)] <- 0   #sm's
    A3[i,(3*n+p+1):(4*n+p)] <- 0   #sp's
  }
  
  b3 <- numeric(n)
  for (i in 1:n) {
    b3[i] <- 0   # for left censoring at 0
  }
  
  A <- rbind(A1,A2,A3)
  
  b <- t(cbind(t(b1),t(b2),t(b3)))

  results3 <- list(c = c, A = A, b = b)
  return(results3)  
}


definelbub <- function(X,y) {
  
  require(pracma)
  #Defines lb,ub for milp.R
  d=50
  #d=40
  #d=20
  #d=10
  #d=5
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  lb1 <-repmat(0,1,n)   #gammas
  lb2 <- repmat(-d,1,p)   #betas
  #lb3 <- repmat(-Inf,1,n)   #phis
  lb3 <- repmat(0,1,n)         #phis, left censoring at zero (0)
  lb4 <- repmat(0,1,n)           #sm's
  lb5 <- repmat(0,1,n)             #sp's
    
  lb <-c(lb1,lb2,lb3,lb4,lb5)
  
  BIGNUM <- Inf
  
  
  ub1 <- repmat(1,1,n)   #gammas
  ub2 <- repmat(d,1,p)     #betas
  ub3 <- repmat(BIGNUM,1,n)  #phis
  ub4 <- repmat(BIGNUM,1,n)    #sm's
  ub5 <- repmat(BIGNUM,1,n)      #sp's

  ub <- c(ub1,ub2,ub3,ub4,ub5)
  

  Aeq <- matrix(0,n,4*n+p)
  for (i in 1:n) {
    Aeq[i,1:n] <- 0   #gammas
    Aeq[i,(n+1):(n+p)] <- 0   #betas
    Aeq[i,(n+p+1):(2*n+p)] <- -1*(1:n==i)  #phis
    Aeq[i,(2*n+p+1):(3*n+p)] <- 1*(1:n==i)   #sm's
    Aeq[i,(3*n+p+1):(4*n+p)] <- -1*(1:n==i)    #sp's
  }
  
  beq <- numeric(n)
  
  for (i in 1:n) {
    beq[i] <- -y[i]
  }
  
  best <- Inf
  
  results4 <- list(lb = lb, ub = ub, Aeq = Aeq, beq = beq, n = n, p = p, best = best)
  return(results4)  
}


milp_cplex <- function(c,A,b,Aeq,beq,lb,ub,n, p, best) {
  
  # Solves a mixed integer lp using cplex 20.1
  # c: is objective function coefficients A: is constraint matrix b: is constraint vector
  # Aeq: is constraint matrix for equality constraints
  # beq: is constraint vector for equality constraints
  # lb: lower bound ub: upper bound n: number of 0-1 variables
  # p: number of predictor variables
  # best: is best solution so far
  # Note this uses the Rcplex Interface documented at
  # https://cran.r-project.org/web/packages/Rcplex/index.html
  # Also, it assumes first n variables must be integer.
  # The MIP equations of the maximum score estimator are available at
  # Bilias, Y., Florios, K., & Skouras, S. (2019). Exact computation of Censored
  # Least Absolute Deviations estimator. Journal of Econometrics, 212(2), 584-606.
  # Writen by Kostas Florios, December 26, 2020
  #
  # C:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex\bin\x64_win64
  #
  # C:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex\include\ilcplex
  #
  
  require(pracma)
  vtype <- c(rep("B",n), rep("C",p), rep("C",n), rep("C",n), rep("C",n))
  
  cvec <- c
  Amat <- rbind(A,Aeq)
  bvec <- c(b,beq)
  Qmat <- NULL
  lb <- lb
  ub <- ub
  control <- list(trace=1,mipemphasis=3,tilim=100)
  #control <- list(tilim=1000)
  #control <- list(tilim=100,mipemphasis=3)
  objsense <- "min"
  sense <- c(rep("L",3*n),rep("E",n))
  
  
  ptm <- proc.time() # Start the clock!
  out <- Rcplex(cvec, Amat, bvec, Qmat = NULL,
                lb = lb, ub = ub, control = control,
                objsense = c("min"), sense = sense, vtype = vtype, n = 1)
  
  elapsed_time <- proc.time() - ptm  # Stop the clock
  cat("Optimization returned status: ",out$status,"\n")
  cat("Objective Value: ",out$obj,"\n")
  cat('(Wall clock) Time elapsed (s): ',elapsed_time[3],"\n")
  cat('Decision variables: show only the betas\n')
  
  x <- out$xopt
  x <- x[(n+1):(n+p)]
  deviation <- out$obj
  feasible <- out$status
  time <- elapsed_time
  results5 <- list(x = x, deviation = deviation, feasible = feasible, time = time)
  
  return(results5)  
}


denormalizeEstimates <- function(estimatesNorm,mu,sigma) {
  
  #denormalized estimatesNorm obtained by Cplex MIP to estimatesRaw, which
  #are meaningful to the user
  
  #quick and dirty implementation, based on GAMS and Fortran Analogues
  p <- length(estimatesNorm)
  betaNorm <- estimatesNorm
  
  betaRaw <- numeric(p)
  betaHelp <- numeric(p)
  
  for (j in 1:p) {
    if (sigma[j] !=0 ) {
      betaHelp[j] <- betaNorm[j] / sigma[j]
    }
    if (sigma[j] ==0 ) {
      for (jj in 1:p) {
        if (sigma[jj] != 0) {
          betaHelp[j] <- betaHelp[j] - betaNorm[jj]*mu[jj]/sigma[jj]
        }
        else {
          jj0 <- jj
        }
      }
      betaHelp[j] <- betaHelp[j] + betaNorm[jj0]
    }
  }
  
  for (j in 1:p) {
    #betaRaw[j] <- betaHelp[j] / betaHelp[1]
    betaRaw[j] <- betaHelp[j]
  }
  
  estimatesRaw <- betaRaw
  results6 <- list(estimatesRaw = estimatesRaw)
  return(results6)  
}


