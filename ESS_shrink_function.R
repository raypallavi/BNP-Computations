#####################################################################
########## Function for MCMC samples using ESS and WC algo ##########
#####################################################################

library(MASS); library(FastGP)

### Function for drawing posterior samples using ESS:

########### For monotone functions estimation ##############
mon.inc.shrink=function(y,x,nu,N,l,eta,mcmc,brn,thin,tau.in,sig.in,xi0.in,xi.in,lam.in,
                    xi0.fix,xi.fix,lam.fix,tau.fix,sig.fix,sseed,verbose,return.plot,return.traplot){
  # y:Response variable; x: vector to form design matrix \Psi (n X N+1)
  # N: (N+1) is the number of knots
  # eta:parameter of the spproximation function of the indicator functions
  # mcmc, brn, thin : mcmc samples, burning and thinning for MCMC
  # tau.in, sig.in, xi0.in, xi.in, lam.in : initial values (supplied by user)
  # xi0.fix,tau.fix,sig.fix,xi.fix : if fixed values of the parameters are to use
  # verbose : logical; if TRUE, prints currnet status; default is TRUE
  # return.plot : logical; if true a plot of estimate with 95% CI is returned; default is TRUE
  # return.traplot : logical; if true traceplots are returned; default is TRUE
  
  # OUTPUT: Posterior samples on xi,lam,xi0,tau,sig and fhat with posterior mean, 95% CI of fhat
  
  if(!missing(lam.fix) && length(unique(lam.fix))==1 && unique(lam.fix)==1)
    return(mon.inc.ESS(y,x,nu,N,l,eta,mcmc,brn,thin,tau.in,sig.in,xi0.in,xi.in,xi0.fix,
                xi.fix,tau.fix,sig.fix,sseed,verbose,return.plot,return.traplot))
  
  if(length(y)!=length(x))
    stop("y and x should be of same length!")
  n=length(y)
  if(missing(N))
    N = ceiling(n/2) - 1
  int_length=1/N
  my_knots=seq(0,1,by=int_length)
  X=des.mat1(x,my_knots,int_length)
  
  if(missing(nu))
    stop("nu needs to be supplied")
  if(nu==0)
    stop("nu cannot be zero")
  if(!missing(l)){
    if(l==0)
      stop("l cannot be zero")
  }
  if(missing(l))
    l=l_est(nu,c(my_knots[1],my_knots[length(my_knots)]),0.05)
  
  # prior covariance K:
  K = covmat(my_knots,nu,l)
  # prior precision:
  K_inv = tinv(K)
  
  if(missing(return.plot))
    return.plot=TRUE
  if(missing(return.traplot))
    return.traplot=TRUE
  
  if(!missing(sseed))
    set.seed(sseed)
  if(missing(sseed))
    set.seed(Sys.Date())
  
  if(missing(verbose))
    verbose=TRUE
  
  if(missing(eta))
    eta=50
  if(missing(mcmc))
    mcmc=5000
  if(missing(brn))
    brn=1000
  if(missing(thin))
    thin=1
  em=mcmc+brn
  ef=mcmc/thin
  
  if(!missing(tau.fix))
    tau.in=tau.fix
  if(!missing(sig.fix))
    sig.in=sig.fix
  if(!missing(xi0.fix))
    xi0.in=xi0.fix
  if(!missing(xi.fix))
    xi.in=xi.fix
  if(!missing(lam.fix))
    lam.in=lam.fix
  
  if(missing(tau.fix) && missing(tau.in))
    tau.in=1
  if(missing(sig.fix) && missing(sig.in))
    sig.in=1
  if(missing(xi0.fix) && missing(xi0.in))
    xi0.in=0
  if(missing(xi.fix) && missing(xi.in))
    xi.in=pmax(0,mvrnorm(1,rep(0,N+1),K))
  if(missing(lam.fix) && missing(lam.in))
    lam.in=rep(1,N+1)
  
  tau=tau.in
  sig=sig.in
  xi0=xi0.in
  xi_in=xi.in
  lam_in=lam.in
  
  gam_in=rep(1,N+1)
  l0=rep(0,ncol(X))
  l1=rep(1,ncol(X))
  
  xi_sam=matrix(0,N+1,ef)
  lam_sam=matrix(0,N+1,ef)
  gam_sam=matrix(0,N+1,ef)
  xi0_sam=rep(0,ef)
  tau_sam=rep(0,ef)
  sig_sam=rep(0,ef)
  fhat_sam=matrix(0,n,ef)
  
  if(verbose)
    print("MCMC sample draws:")
  
  ptm=proc.time()
  for(i in 1:em){
    # sampling Xi:
    y_tilde = y - xi0
    if(missing(xi.fix)){
      nu.xi = as.vector(samp.WC(my_knots,nu,l,tau))
      xi_out = ESS(xi_in,nu.xi,y_tilde,t(t(X)*lam_in),sig,eta)
      xi_out = sapply(xi_out,function(z) return(max(0,z)))
    }else{
      xi_in = xi_out
    }
    
    # sampling lambda:
    if(missing(lam.fix)){
      nu.lam = rnorm(l1,l0,sd=1/sqrt(gam_in))
      lam_out = ESS(lam_in,nu.lam,y_tilde,t(t(X)*xi_out),sig,eta)
      lam_out = sapply(lam_out,function(z) return(max(0,z)))
    }else{
      lam_in = lam_out
    }
    
    # sampling xi_0:
    Xxi = X%*%(xi_out*lam_out)
    y_star = y - Xxi
    if(missing(xi0.fix)){
      xi0 = rnorm(1,mean(y_star),sqrt(sig/n))
    }
    
    # sampling gamma_j's:
    if(missing(lam.fix))
      gam_out = rexp(l1,rate = 0.5*(1+lam_out^2))
    
    # sampling \sigma^2:
    y0 = y_star - xi0
    if(missing(sig.fix))
      sig = 1/rgamma(1,shape = n/2,rate = sum(y0^2)/2)
    
    # sampling \tau^2 from horseshoe code:
    if(missing(tau.fix)){
      tempt = (t(xi_out)%*%K_inv%*%xi_out)/2 
      et = 1/tau^2
      utau = runif(1,0,1/(1+et))
      ubt = min(1,(1-utau)/utau)
      Fubt = pgamma(ubt,((N+1)+1)/2,scale=1/tempt) 
      Fubt = max(Fubt,1e-10) # for numerical stability
      ut = runif(1,0,Fubt)
      et = qgamma(ut,((N+1)+1)/2,scale=1/tempt) 
      tau = 1/sqrt(et)
    }
    
    # storing MCMC samples:
    if(i > brn && i%%thin == 0){
      xi_sam[,(i-brn)/thin]=xi_out
      lam_sam[,(i-brn)/thin]=lam_out
      gam_sam[,(i-brn)/thin]=gam_out
      xi0_sam[(i-brn)/thin]=xi0
      sig_sam[(i-brn)/thin]=sig
      tau_sam[(i-brn)/thin]=tau
      fhat_sam[,(i-brn)/thin]=xi0 + Xxi
    }
    
    if(i%%1000==0 && verbose){
      print(i)
    }
    
    # renewing the intial value:
    xi_in = xi_out; lam_in = lam_out; gam_in=gam_out
    
  }; tm=proc.time()-ptm
  
  if(return.traplot){
    library(mcmcplots)
    mx=min(N+1,25)
    xi_mat=matrix(xi_sam[1:mx,],nrow = mx,
                  dimnames = list(c(paste("xi[",1:mx,"]",sep = "")),NULL))
    traplot(t(xi_mat))
    
    lam_mat=matrix(lam_sam[1:mx,],nrow = mx,
                  dimnames = list(c(paste("lambda[",1:mx,"]",sep = "")),NULL))
    traplot(t(lam_mat))
    
    par(mfrow=c(1,1))
    traplot(as.mcmc(xi0_sam),main=expression(paste(xi[0])))
    traplot(as.mcmc(sig_sam),main=expression(paste(sigma^2)))
    traplot(as.mcmc(tau_sam),main=expression(paste(tau^2)))
  }
  
  fmean=rowMeans(fhat_sam)
  qnt=apply(fhat_sam,1,function(x) quantile(x,c(0.025,0.975),na.rm = TRUE))
  f_low=qnt[1,]
  f_upp=qnt[2,]
  ub=max(f_low,f_upp,fmean)
  lb=min(f_low,f_upp,fmean)
  
  if(return.plot){
    par(mfrow=c(1,1))
    plot(x,fmean,type="l",lwd=2,main="Black: Estimate, Blue: 95% CI",ylim=c(lb,ub))
    lines(x,f_low,lwd=2,lty=2,col="blue")
    lines(x,f_upp,lwd=2,lty=2,col="blue")
  }
  
  return(list("time"=tm,"xi_sam"=xi_sam,"lam_sam"=lam_sam,"xi0_sam"=xi0_sam,"sig_sam"=sig_sam,
              "tau_sam"=tau_sam,"fhat_sam"=fhat_sam,"fmean"=fmean,"f_low"=f_low,"f_upp"=f_upp))
}




#END