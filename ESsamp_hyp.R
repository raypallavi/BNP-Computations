### Code that implements simulations discussed in the paper:

# source codes and libraries:
require(Rcpp); library(ggplot2)
source("all_functions.R")
sourceCpp("inv_chol.cpp")
source("ESS_hyp_functions.R")

## A monotone function (as in Maatouk & Bay)
f1=function(x){
  return(log(20*x + 1))
}

f2=function(x){
  return(5 * (x-0.5)^2)
}

#####################################################################
########## Monotone (shape-constraint) function estimation ##########
#####################################################################

#### generating the data ####
n=500
x=sort(runif(n,min=0,max=1))
y.true=sapply(x,FUN = f1)
### finding noise variance:
signal = y.true
noise = rnorm(n,sd = 0.7)
### Y-data
y=signal+noise
print(sd(signal)/sd(noise))

ntr=300; nte=n-ntr
ind=sort(sample(1:n,ntr))
xtrain=x[ind]
ytrain=y[ind]
xtest=x[-ind]
ytest=y[-ind]
yte.tr=y.true[-ind]

# knots and design matrix:
N=ceiling(ntr/2)-1
int_length=1/N
my_knots=seq(0,1,by=int_length)

Xtr=matrix(0,ntr,N+1)                   # design matrix \psi(n X N+1)
for(j in 1:ntr){
  Xtr[j,1:(N+1)]=psi_j(xtrain[j],my_knots,int_length)
}
#XtrTXtr=t(Xtr)%*%Xtr

Xte=matrix(0,nte,N+1)                   # design matrix \psi(n X N+1)
for(j in 1:nte){
  Xte[j,1:(N+1)]=psi_j(xtest[j],my_knots,int_length)
}

# Posterior samples:
post.sam = mon.inc.ESS.hyp(ytrain,xtrain,eta=50,mcmc=5000,brn=2000,thin=1,
                           sseed=Sys.Date(),return.plot = FALSE,return.traplot = FALSE)

### Prediction:
xi_sam=post.sam$xi_sam
ef=ncol(xi_sam)
xi0_sam=post.sam$xi0_sam
fhat_sam=t(xi0_sam + t(Xte %*% xi_sam))
fmean=rowMeans(fhat_sam)
fsort=fhat_sam
for(i in 1:nrow(fhat_sam)){
  fsort[i,]=sort(fhat_sam[i,])
}
f_low=fsort[,ef*0.025]
f_upp=fsort[,ef*0.975]

### using gg-plot:
pl = data.frame(xtest,ytest,fmean)
pl$menlelb = f_low
pl$menleub = f_upp
p<-ggplot(pl, aes(x=xtest,y=ytest)) +
  geom_ribbon(aes(ymin=menlelb,ymax=menleub),alpha=0.5,show.legend=TRUE) +
  geom_line(aes(y=fmean), colour="blue",size=1) +
  geom_line(aes(y=yte.tr), colour="red", linetype=2,size=1) + 
  geom_point(colour="green3",size=1.5)+labs(x = "", y="")
df=data.frame(xtest=xtrain,ytest=ytrain)
p+geom_point(data = df,colour="green3",size=1,shape=24,fill="green3",stroke=1)

#####################################################################
########### Convex (shape-constraint) function estimation ###########
#####################################################################

#### generating the data ####
n=500
x=sort(runif(n,max=1))
y.true=sapply(x,FUN = f2)
### finding noise variance:
signal = y.true
noise = rnorm(n,sd = 0.7)
### Y-data
y=signal+noise
print(sd(signal)/sd(noise))
ntr=300; nte=n-ntr
ind=sort(sample(1:n,ntr))
xtrain=x[ind]
ytrain=y[ind]
xtest=x[-ind]
ytest=y[-ind]
yte.tr=y.true[-ind]

# knots and design matrix:
N=ceiling(ntr/2)-1
int_length=1/N
my_knots=seq(0,1,by=int_length)

Xtr=matrix(0,ntr,N+1)                   # design matrix \phi(n X N+1)
for(j in 1:ntr){
  Xtr[j,1:(N+1)]=phi_j(xtrain[j],my_knots,int_length)
}
#XtrTXtr=t(Xtr)%*%Xtr

Xte=matrix(0,nte,N+1)                   # design matrix \phi(n X N+1)
for(j in 1:nte){
  Xte[j,1:(N+1)]=phi_j(xtest[j],my_knots,int_length)
}


# Posterior samples:
post.sam=con.ESS.hyp(ytrain,xtrain,eta=50,mcmc=5000,brn=2000,thin=1,
                     sseed=Sys.Date(),return.plot = FALSE,return.traplot = FALSE)

### Prediction:
xi_sam=post.sam$xi_sam
ef=ncol(xi_sam)
xi0_sam=post.sam$xi0_sam
xi1_sam=post.sam$xi1_sam

z=mapply(FUN=xix,a=xi1_sam,MoreArgs = list(b=xtest))
fhat_sam=t(xi0_sam + t(z) + t(Xte %*% xi_sam))

fmean=rowMeans(fhat_sam)
fsort=fhat_sam
for(i in 1:nrow(fhat_sam)){
  fsort[i,]=sort(fhat_sam[i,])
}
f_low=fsort[,ef*0.025]
f_upp=fsort[,ef*0.975]

### using gg-plot:
pl = data.frame(xtest, fmean)
pl$menlelb = f_low
pl$menleub = f_upp
p<-ggplot(pl, aes(x=xtest,y=ytest)) +
  geom_ribbon(aes(ymin=menlelb,ymax=menleub),alpha=0.5,show.legend=TRUE) +
  geom_line(aes(y=fmean), colour="blue",size=1) +
  geom_line(aes(y=yte.tr), colour="red", linetype=2,size=1) + 
  geom_point(colour="green3",size=1.5)+labs(x = "", y="")
df=data.frame(xtest=xtrain,ytest=ytrain)
p+geom_point(data = df,colour="green3",size=1,shape=24,fill="green3",stroke=1)



# END OF SIMULATION CODE