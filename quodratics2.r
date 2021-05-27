library(expm)
#library(pracma)

#Covariance:
##[[1            0.999999999]
## [0.999999999            1]]

rho <- 1-1/10^7
U<-function(q) .5 *(q[1]^2 +q[2]^2-2*q[1]*q[2]*rho)/(1-rho^2)
dU<-function(q) c(.5*(2*q[1]-2*q[2]*rho)/(1-rho^2),.5*(2*q[2]-2*q[1]*rho)/(1-rho^2))
#dU<-function(q) grad(U,q)
mat<-matrix(c(1/(1-rho^2),-rho/(1-rho^2),-rho/(1-rho^2),1/(1-rho^2)),nrow=2,byrow=T)
ddU<-function(q) mat

##ddU<-function(q) hessian(U,q,.Machine$double.eps^(1/4))

source("./Sampler1.r")
BURNIN <- 500
EPISODE <- 1000
Dim <- 2
CHAINS <- 50
tic()
QS<-array(hmc(U,dU,ddU,Dim,BURNIN,EPISODE,T,T, NA),dim=c(Dim,CHAINS*(EPISODE-BURNIN)))
toc(echo=TRUE)
plot(t(QS))
QS1 <- sqrtm(mat) %*% QS
plot(t(QS1))
