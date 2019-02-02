# Source files and libraries:
require(Rcpp)
source("all_functions.R")
source("ESS_functions.R")
sourceCpp("inv_chol.cpp")
source("ESS_hyp_functions.R")

# Call data .csv file
dat=read.csv("cars-mbart.csv",header = TRUE,sep = ",")
x1=-dat$mileage      #since the data is on negative mileage
xmax=max(x1)        #scaled to lie within [0,1]
x1=x1/xmax; x=sort(x1)
y1=dat$price; y=log(y1[order(x1)]) #log-transform of sale-price
plot(x,y,type = "p",pch=20)

###### With fixed hyperparameters ######
## calculating MSPE for one split of the 10 splits
n=length(y); ntr=900; nte=n-ntr
ind=sort(sample(1:n,ntr))
xtrain=x[ind]; ytrain=y[ind]
xtest=x[-ind]; ytest=y[-ind]

## Using monotone decreasing criteria and without hyperparameter updates:
post.sam=mon.dec.ESS(ytrain,xtrain,nu=0.5,eta=50,mcmc=5000,brn=1000,thin=1,
                     sseed=Sys.Date(),return.plot=FALSE,return.traplot=FALSE)
# Prediction using test data:
xi_sam=post.sam$xi_sam
xi0_sam=post.sam$xi0_sam
N=nrow(xi_sam)-1; int_length=1/N; my_knots=seq(0,1,by=int_length)
Xte=matrix(0,nte,N+1)                   # design matrix \psi(n X N+1)
for(j in 1:nte){
  Xte[j,1:(N+1)]=psi_j(xtest[j],my_knots,int_length)
}
fhat_sam=t(xi0_sam + t(Xte %*% xi_sam))
fmean=rowMeans(fhat_sam)
(mspe=mean((fmean - ytest)^2))

###### Without fixing hyperparameters ######
## calculating MSPE for one split of the 10 splits
n=length(y); ntr=900; nte=n-ntr
ind=sort(sample(1:n,ntr))
xtrain=x[ind]; ytrain=y[ind]
xtest=x[-ind]; ytest=y[-ind]

## Using monotone decreasing criteria and with hyperparameter updates:
post.sam=mon.dec.ESS.hyp(ytrain,xtrain,eta=50,sseed=1234,mcmc=5000,brn=1000,thin=1,
                     return.plot=FALSE,return.traplot=FALSE)
# Prediction on test data:
xi_sam=post.sam$xi_sam
xi0_sam=post.sam$xi0_sam
N=nrow(xi_sam)-1; int_length=1/N; my_knots=seq(0,1,by=int_length)
Xte=matrix(0,nte,N+1)                   # design matrix \psi(n X N+1)
for(j in 1:nte){
  Xte[j,1:(N+1)]=psi_j(xtest[j],my_knots,int_length)
}
fhat_sam=t(xi0_sam + t(Xte %*% xi_sam))
fmean=rowMeans(fhat_sam)
(mspe=mean((fmean - ytest)^2))

# END