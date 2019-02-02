### All required functions for using ESS
### Functions related to Wood and Chan algorithm of drawing samples
### MH for sampling from \nu and \ell
### Covariance matrix and design matrix (using basis function) are also defined
### And all related and dependant functions are here

### Required libraries:
library(fields)
library(FastGP)

### Define design matrix ###
### The basis functions ###

# For monotone function estimation:
psi_j=function(x,my_knot,delta_N)
{
  N=length(my_knot)
  k=rep(0,N)
  i=max(which(my_knot<=x))
  if(i==1)
  {
    k[1]=x-0.5*(x^2)/delta_N
    k[2]=x-my_knot[2]*x/delta_N+0.5*x^2/delta_N
  }
  if(i==2)
  {
    k[1]=delta_N/2
    k[2]=delta_N/2+(x-my_knot[2])*(1+my_knot[2]/delta_N)-0.5*(x^2-my_knot[2]^2)/delta_N
    k[3]=(x-my_knot[2])*(1-my_knot[3]/delta_N)+0.5*(x^2-my_knot[2]^2)/delta_N
  }
  if(i==N)
  {
    k[1]=delta_N/2
    k[2:(N-1)]=delta_N
    k[N]=delta_N/2
  }
  if(i!=1 && i!=2 && i!=N)
  {
    k[1]=delta_N/2
    k[2:(i-1)]=delta_N
    k[i]=delta_N/2+(x-my_knot[i])*(1+my_knot[i]/delta_N)-0.5*(x^2-my_knot[i]^2)/delta_N
    k[i+1]=(x-my_knot[i])*(1-my_knot[i+1]/delta_N)+0.5*(x^2-my_knot[i]^2)/delta_N
  }
  return(k)
}

# For convex function estimation:
phi_j=function(x,my_knot,delta_N)
{
  N=length(my_knot)
  k=rep(0,N)
  if(x>=my_knot[1] && x<my_knot[2])
  {
    k[1]=(x^2/2)-0.5*(x^3/3)/delta_N
    k[2]=0.5*(x^3/3)/delta_N
  }
  if(x>=my_knot[2] && x<my_knot[3])
  {
    k[1]=(my_knot[2]^2/2)-0.5*(my_knot[2]^3/3)/delta_N+0.5*delta_N*(x-my_knot[2])
    k[2]=my_knot[2]^3/(3*delta_N)+0.5*delta_N*(x-my_knot[2])+(x^2-my_knot[2]^2)-1.5*my_knot[2]*(x-my_knot[2])-0.5*(x^3/3)/delta_N
    k[3]=(1-my_knot[3]/delta_N)*(x^2/2-x*my_knot[2]+my_knot[2]^2/2)+0.5*(x^3/3-x*my_knot[2]^2+2*my_knot[2]^3/3)/delta_N
  }
  if(x>=my_knot[3] && x<my_knot[N])
  {
    k[1]=(my_knot[2]^2/2)-0.5*(my_knot[2]^3/3)/delta_N+0.5*delta_N*(x-my_knot[2])
    k[2]=my_knot[2]^3/(3*delta_N)+0.5*delta_N^2+my_knot[3]^2-2.5*my_knot[2]^2-0.5*(my_knot[3]^3/3)/delta_N+delta_N*(x-my_knot[3])
    k[3]=(1-my_knot[3]/delta_N)*(my_knot[3]^2/2-my_knot[3]*my_knot[2]+my_knot[2]^2/2)+0.5*(my_knot[3]^3/3-my_knot[3]*my_knot[2]^2+2*my_knot[2]^3/3)/delta_N+delta_N*(x-my_knot[3])
    
    if(N>4){
      for(j in 4:(N-1)){
        if(x<my_knot[j-1])
          k[j]=0
        if(x>=my_knot[j-1] && x<my_knot[j])
          k[j]=(1-my_knot[j]/delta_N)*(0.5*(x^2-my_knot[j-1]^2)-my_knot[j-1]*(x-my_knot[j-1]))+0.5*((x^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(x-my_knot[j-1]))/delta_N
        if(x>=my_knot[j] && x<my_knot[j+1])
          k[j]=(1-my_knot[j]/delta_N)*(0.5*(my_knot[j]^2-my_knot[j-1]^2)-my_knot[j-1]*(my_knot[j]-my_knot[j-1]))+0.5*((my_knot[j]^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(my_knot[j]-my_knot[j-1]))/delta_N+0.5*delta_N*(x-my_knot[j])+
            (1+my_knot[j]/delta_N)*(0.5*(x^2-my_knot[j]^2)-my_knot[j]*(x-my_knot[j]))-0.5*((x^3/3-my_knot[j]^3/3)-my_knot[j]^2*(x-my_knot[j]))/delta_N
        if(x>=my_knot[j+1])
          k[j]=(1-my_knot[j]/delta_N)*(0.5*(my_knot[j]^2-my_knot[j-1]^2)-my_knot[j-1]*(my_knot[j]-my_knot[j-1]))+0.5*((my_knot[j]^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(my_knot[j]-my_knot[j-1]))/delta_N+0.5*delta_N*(my_knot[j+1]-my_knot[j])+
            (1+my_knot[j]/delta_N)*(0.5*(my_knot[j+1]^2-my_knot[j]^2)-my_knot[j]*(my_knot[j+1]-my_knot[j]))-0.5*((my_knot[j+1]^3/3-my_knot[j]^3/3)-my_knot[j]^2*(my_knot[j+1]-my_knot[j]))/delta_N+delta_N*(x-my_knot[j+1])
      }
    }
    
    if(x>=my_knot[N-1] && x<my_knot[N])
      k[N]=(1-my_knot[N]/delta_N)*(0.5*(x^2-my_knot[N-1]^2)-my_knot[N-1]*(x-my_knot[N-1]))+0.5*((x^3/3-my_knot[N-1]^3/3)-my_knot[N-1]^2*(x-my_knot[N-1]))/delta_N
  }
  if(x>=my_knot[N])
  {
    k[1]=(my_knot[2]^2/2)-0.5*(my_knot[2]^3/3)/delta_N+0.5*delta_N*(x-my_knot[2])
    k[2]=my_knot[2]^3/(3*delta_N)+0.5*delta_N^2+my_knot[3]^2-2.5*my_knot[2]^2-0.5*(my_knot[3]^3/3)/delta_N+delta_N*(x-my_knot[3])
    k[3]=(1-my_knot[3]/delta_N)*(my_knot[3]^2/2-my_knot[3]*my_knot[2]+my_knot[2]^2/2)+0.5*(my_knot[3]^3/3-my_knot[3]*my_knot[2]^2+2*my_knot[2]^3/3)/delta_N+delta_N*(x-my_knot[3])
    for(j in 4:(N-1)){
      k[j]=(1-my_knot[j]/delta_N)*(0.5*(my_knot[j]^2-my_knot[j-1]^2)-my_knot[j-1]*(my_knot[j]-my_knot[j-1]))+0.5*((my_knot[j]^3/3-my_knot[j-1]^3/3)-my_knot[j-1]^2*(my_knot[j]-my_knot[j-1]))/delta_N+0.5*delta_N*(my_knot[j+1]-my_knot[j])+
        (1+my_knot[j]/delta_N)*(0.5*(my_knot[j+1]^2-my_knot[j]^2)-my_knot[j]*(my_knot[j+1]-my_knot[j]))-0.5*((my_knot[j+1]^3/3-my_knot[j]^3/3)-my_knot[j]^2*(my_knot[j+1]-my_knot[j]))/delta_N+delta_N*(x-my_knot[j+1])
    }
    k[N]=(1-my_knot[N]/delta_N)*(0.5*(my_knot[N]^2-my_knot[N-1]^2)-my_knot[N-1]*(my_knot[N]-my_knot[N-1]))+0.5*((my_knot[N]^3/3-my_knot[N-1]^3/3)-my_knot[N-1]^2*(my_knot[N]-my_knot[N-1]))/delta_N+0.5*delta_N*(x-my_knot[N])
  }
  return(k)
}

### Function to form design matrix:
des.mat1=function(x,my_knot,delta_N){
  # Function to form basis matrix for monotone constraint
  n=length(x)
  N=length(my_knot)-1
  # design matrix \Psi(n X N+1)
  X=matrix(0,n,N+1)
  for(l in 1:n){
    X[l,1:(N+1)]=psi_j(x[l],my_knot,delta_N)
  }
  return(X)
}

des.mat2=function(x,my_knot,delta_N){
  # Function to form basis matrix for convex constraint
  n=length(x)
  N=length(my_knot)-1
  # design matrix \Phi(n X N+1)
  X=matrix(0,n,N+1)
  for(l in 1:n){
    X[l,1:(N+1)]=phi_j(x[l],my_knot,delta_N)
  }
  return(X)
}

### Function to compute \xi_1 * x :
xix=function(a,b){
  return(a*b)
}

#Given a \nu (smoothness parameter of matern kernel) finding a value of 
# l (length-scale parameter) such that the correlation between the 
# maximum seperation is some small value, say 0.05

#Matern kernel with smoothness nu and length-scale l:
MK = function(x, y ,l, nu){
  ifelse(abs(x-y)>0, (sqrt(2*nu)*abs(x-y)/l)^nu/(2^(nu-1)*gamma(nu))*besselK(x=abs(x-y)*sqrt(2*nu)/l, nu=nu), 1.0)
}

# function for uniroot:
fl=function(l,para){ 
  #para[1]=x, para[2]=y and para[3]=nu of MK : Matern kernel function;
  #para[4]=pre-specified value of the correlation
  a=MK(para[1],para[2],l,para[3])
  return(a-para[4])
}

# function for estimating l:
l_est=function(nu,range,val){
  # nu : smoothness; range : c(min, max) of the range of variable
  # val : pre-specified value of the correlation between the maximum seperation
  para=c(range[1],range[2],nu,val)
  rl=uniroot(f=fl,interval = c(0.000001,100000),para)
  return(rl$root)
}

# Covariance matrix
covmat=function(knot,nu,l){
  return(MK(rdist(knot),0,l,nu))
}

# Order of the circulant matrix:
# minimum value of g and m so that G can be embedded into C
min_g=function(knot){
  N=length(knot)
  g=ceiling(log(2*N,2))   #m=2^g and m>=2(n-1) : Wood & Chan notation; 
  #since we are going upto n and not stopping at (n-1), the condition is modified!
  return("g" = g)
}

# forming the circulant matrix:
circulant=function(x){
  n = length(x)
  mat = matrix(0, n, n)
  for (j in 1:n) {
    mat[j, ] <- c(x[-(1:(n+1-j))], x[1:(n+1-j)])
  }
  return(mat)
}

# Function for forming the vector of circulant matrix:
circ_vec=function(knot,g,nu,l,tausq){
  delta_N=1/(length(knot)-1)
  m=2**g
  cj=integer()
  for(j in 1:m){
    if(j<=(m/2))
      cj[j]=(j-1)*delta_N
    else
      cj[j]=(m-(j-1))*delta_N
  }
  x=(tausq*MK(cj,0,l,nu))
  return(x)
}

# Function for finding a g such that C is nnd:
eig.eval=function(knot,g,nu,l,tausq){
  vec=circ_vec(knot,g,nu,l,tausq)
  C=circulant(vec)
  ev=min(eigen(C)$values)
  return(list("vec" = vec, "min.eig.val" = ev))
}

# Function for finding a g such that C is nnd:
# without forming the circulant matrix and without computing eigen values:
C.eval=function(knot,g,nu,l,tausq){
  vec=circ_vec(knot,g,nu,l,tausq)
  val=fft(vec) # eigenvalues will be real as the circulant matrix formed by the 
               # vector is by construction is symmetric!
  ev=min(Re(val))
  return(list("vec" = vec, "min.eig.val" = ev))
}


nnd_C=function(knot,g,nu,l,tausq){
  C.vec=C.eval(knot,g,nu,l,tausq)$vec
  eval=C.eval(knot,g,nu,l,tausq)$min.eig.val
  if(eval>0)
    return(list("cj" = C.vec,"g" = g))
  else{
    g=g+1
    nnd_C(knot,g,nu,l,tausq)
  }
}

# computing the eigen values of C using FFT:
eigval=function(knot,nu,l,tausq){
  g=min_g(knot)
  c.j=nnd_C(knot,g,nu,l,tausq)$cj
  lambda=Re(fft(c.j))
  if(min(lambda)>0)
    return(lambda)
  else
    stop("nnd condition is NOT satisfied!!!")
}


#################################################################
########## Samples drawn using Wood and Chan Algorithm ##########
#################################################################
samp.WC=function(knot,nu,l,tausq){
  N=length(knot)
  lambda=eigval(knot,nu,l,tausq)
  m=length(lambda)
  samp.vec=rep(0,N)
  a=rep(0,m)
  a[1]=sqrt(lambda[1])*rnorm(1)/sqrt(m)
  a[(m/2)+1]=sqrt(lambda[(m/2)+1])*rnorm(1)/sqrt(m)
  i=sqrt(as.complex(-1))
  for(j in 2:(m/2)){
    uj=rnorm(1); vj=rnorm(1)
    a[j]=(sqrt(lambda[j])*(uj + i*vj))/(sqrt(2*m))
    a[m+2-j]=(sqrt(lambda[j])*(uj - i*vj))/(sqrt(2*m))
  }
  samp=fft(a)
  samp.vec=Re(samp[1:N])
  return(samp.vec)
}

#############################################
########## Functions for using ESS ##########
#############################################
ESS = function(beta,nu_ess,y,X,sigsq,eta){
  thetamin = 0; 
  thetamax = 2*pi;
  
  u = runif(1)
  logy = loglik(y,X,sigsq,eta,beta) + log(u); 
  
  theta = runif(1,thetamin,thetamax); 
  thetamin = theta - 2*pi; 
  thetamax = theta;
  betaprime = beta*cos(theta) + nu_ess*sin(theta);
  
  while(loglik(y,X,sigsq,eta,betaprime) <= logy){
    if(theta < 0)
      thetamin = theta
    else
      thetamax = theta
    theta = runif(1,thetamin,thetamax)
    betaprime = beta*cos(theta) + nu_ess*sin(theta)
  }
  return(betaprime)       
}

ESS.dec = function(beta,nu_ess,y,X,sigsq,eta){
  thetamin = 0; 
  thetamax = 2*pi;
  
  u = runif(1)
  logy = loglik2(y,X,sigsq,eta,beta) + log(u); 
  
  theta = runif(1,thetamin,thetamax); 
  thetamin = theta - 2*pi; 
  thetamax = theta;
  betaprime = beta*cos(theta) + nu_ess*sin(theta);
  
  while(loglik2(y,X,sigsq,eta,betaprime) <= logy){
    if(theta < 0)
      thetamin = theta
    else
      thetamax = theta
    theta = runif(1,thetamin,thetamax)
    betaprime = beta*cos(theta) + nu_ess*sin(theta)
  }
  return(betaprime)       
}

## Defining the loglik function to be used in ESS:
## loglik calculates the log of the likelihood:
loglik=function(y,X,sigsq,eta,beta){
  mu=y-(X%*%beta)
  val=eta*sum(beta)-sum(log(1+exp(eta*beta)))-sum(mu^2)/(2*sigsq)
  return(val)
}

loglik2=function(y,X,sigsq,eta,beta){
  mu=y-(X%*%beta)
  val=-sum(log(1+exp(eta*beta)))-sum(mu^2)/(2*sigsq)
  return(val)
}


## MH algo for \nu of Matern kernel: (NOT USED IN THE MAIN CODE)
nu.MH1 = function(nu.in,l.in,tau.in,xi.in,Kmat,knot,range.nu=c(0.5,1),sd.p=0.05){
  nu.cand = exp(log(nu.in)+rnorm(1,0,sd.p))
  l.cand = l_est(nu.cand,c(0,1),0.05)
  du = dunif(nu.cand,range.nu[1],range.nu[2])
  if(du > 0){
    Kcand = covmat(knot,nu.cand,l.cand)
    Linv = inv_chol(Kmat); Linv.cand = inv_chol(Kcand)
    r = exp(sum(log(diag(Linv.cand)))-sum(log(diag(Linv))))*(nu.cand/nu.in)
    t1 = sum((t(Linv.cand)%*%xi.in)^2); t2 = sum((t(Linv)%*%xi.in)^2)
    r = r*exp(-(t1 - t2)/(2*tau.in))
    alpha = min(r,1)
  }
  else
    alpha = 0
  u = runif(1)
  nu.out = (u < alpha)*nu.cand + (u >= alpha)*nu.in
  l.out = (u < alpha)*l.cand + (u >= alpha)*l.in
  cnt = as.numeric((abs(nu.out - nu.in) > 0))
  return(list("nu" = nu.out,"l" = l.out,"count" = cnt))
}

## MH algo for \nu and \ell of Matern kernel:
nu.MH2 = function(nu.in,l.in,tau.in,xi.in,knot,range.nu=c(0.5,1),range.l=c(0.1,1),sd.nu=0.05,sd.l=0.1){
  Kmat = covmat(knot,nu.in,l.in)
  Linv = inv_chol(Kmat)
  nu.cand = exp(log(nu.in)+rnorm(1,0,sd.nu))
  l.cand = exp(log(l.in)+rnorm(1,0,sd.l))
  dnu = dunif(nu.cand,range.nu[1],range.nu[2])
  dl = dunif(l.cand,range.l[1],range.l[2])
  if(dnu > 0 && dl > 0){
    Kcand = covmat(knot,nu.cand,l.cand)
    Linv.cand = inv_chol(Kcand)
    t1 = sum((t(Linv.cand)%*%xi.in)^2)
    t2 = sum((t(Linv)%*%xi.in)^2)
    r = exp(sum(log(diag(Linv.cand)))-sum(log(diag(Linv)))-((t1 - t2)/(2*tau.in)))*(nu.cand/nu.in)*(l.cand/l.in)
    alpha = min(r,1)
  }
  else{
    alpha = 0
    Linv.cand = Linv
  }
  u = runif(1)
  nu.out = (u < alpha)*nu.cand + (u >= alpha)*nu.in
  l.out = (u < alpha)*l.cand + (u >= alpha)*l.in
  cnt = (u < alpha)*1 + (u >= alpha)*0
  L_inv = (u < alpha)*Linv.cand + (u >= alpha)*Linv
  return(list("nu" = nu.out,"l" = l.out,"count" = cnt,"L_inv"=L_inv))
}


