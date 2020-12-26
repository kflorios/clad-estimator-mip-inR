require(Rcplex)

source("functions.R")

#Main Program: Computes CLAD by defining a MILP and calling milp function

results1 <- readXyw()
X <- results1$X
y <- results1$y
w <- results1$w
results2 <- standardizeX(X)
X  <- results2$X
mu <- results2$mu
sigma <- results2$sigma
results3 <- definecAb(X,y,w)
c <- results3$c
A <- results3$A
b <- results3$b
results4 <- definelbub(X,y)
lb <- results4$lb
ub <- results4$ub
Aeq <- results4$Aeq
beq <- results4$beq
n <- results4$n
p <- results4$p
best <- results4$best
results5 <- milp_cplex(c,A,b,Aeq,beq,lb,ub,n,p,best)
x <- results5$x
deviation <- results5$deviation
feasible <- results5$feasible
time <- results5$time
estimatesNorm=x
value=deviation
status=feasible
runtime=time
quality=status

estimatesRaw=denormalizeEstimates(estimatesNorm,mu,sigma)
estimates=estimatesRaw
