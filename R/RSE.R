#RSE
#rare species estimation package
#####################################


################################################
#general functions
X.to.f = function(X){
  f = factor(X, levels = 0:max(X))
  f = table(f, exclude = 0) ## frequency counts
  dimnames(f) = NULL
  return(f)
} 

f.to.X = function(f){
  X = c()
  k=1
  while(k<=length(f)){
    if(k>1) X = c(X,rep(k,f[k])) else {
      X = rep(k,f[k])
    }
    k = k+1;
  }
  return(X);
} 


################################################
#For abundance-based data
### Chao et al. (2015) Ecology paper: species-Rank abundance distribution
DetAbu <- function(x, zero=FALSE){
  x <- unlist(x)
  # print(x)
  n <- sum(x)  
  f1 <- sum(x==1)
  f2 <- sum(x==2)
  f3 <- sum(x==3)
  if(f2==0){
    f1 <- max(f1 - 1, 0)
    f2 <- 1
  }
  A1 <- f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*max(f2,1)))
  A2 <- f2 / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3*max(f3,1)))^2
  if(zero==FALSE) x <- x[x>0]
  q.solve <- function(q){
    e <- A1 / sum(x/n*exp(-q*x))
    out <- sum((x/n * (1 - e * exp(-q*x)))^2) - sum(choose(x,2)/choose(n,2)) + A2
    abs(out)
  }
  q <- tryCatch(optimize(q.solve, c(0,1))$min, error = function(e) {1})
  e <- A1 / sum(x/n*exp(-q*x))
  o <- x/n * (1 - e * exp(-q*x))
  return(c(e,q))
  # o
}


### individual-based Chao et al. (2013): 
### Generating species counts by the bootstrapping method
boot.abundance.fun = function(S.hat,f,b) {
  
  D = sum(f)
  ind = 1:length(f)
  kmax = length(f)
  d = rep(0,kmax)
  n = sum(ind*f)
  S.hat = ceiling(S.hat)
  p0=0;f0=0;
  
  if(f[2]==0){
    f1 <- max(f[1] - 1, 0)
    f2 <- 1
  } else {
    f1 = f[1]
    f2 = f[2]
  }
  C = 1-f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*max(f2,1)))
  
  d = ind/n*(1-b[1]*exp(-b[2]*ind))
  
  cell.p = rep(d, times=f)
  label.X = 1:length(cell.p)

  if(S.hat>D) {
    f0=S.hat-D;
    p0=(1-C)/f0
    sp.count = sample(c(length(label.X)+(1:f0),label.X), n, 
                      replace = T, prob = c(rep(p0,f0),cell.p))
  } else {
    sp.count = sample(label.X, n, replace = T, prob = c(cell.p))
  }
  sp.count = factor(sp.count, levels = 1:(length(label.X)+f0))

  sp.count =table(sp.count)

  f.count = X.to.f(sp.count)

  return(f.count)
}



##### The Bayesian-weight method
Pred.Fk.BW = function(f,m,b,k.show=3){
  kmax = length(f)
  
  cut.pts = max(10,kmax)
  
  n = sum(f*(1:kmax))
  ## delta k
  d = rep(0,cut.pts) 
  ## estimated fk
  est.f = rep(0,k.show) 
  

  ### relative abundance estimation by Chao et al. (2015)
  ind = 1:kmax
  d = ind/n*(1-b[1]*exp(-b[2]*ind))
  
  for(i in 1:k.show){## i for k
    for(j in 1:min(i,kmax)){## j for q
      # est.f[i] = est.f[i] +
      #   f[j]*choose(m,i)/(choose(m+n,i)-choose(m,i))*
      #   dbinom(i-j,m,d[j]);
      
      tmp = prod(1+n/(m-(0:(i-1))))
      est.f[i] = est.f[i] +
        f[j]*1/(tmp-1)*dbinom(i-j,m,d[j]);
    }
  }
  
  return(est.f)
}

### Chao1 estimator
SpEst.Chao1.abun = function(f) {
  est = sum(f)
  n = sum(f * (1:length(f)))
  if (length(f) > 1) {
    if (f[2] > 0) {
      A = f[1]/f[2]
    } else {
      if (f[1] > 0) {
        A = f[1] - 1
      } else {
        A = 0
      }
    }
  } else A = f[1] - 1
  est = est + (n - 1)/n * A * f[1]/2
  return(est)
}



##### the naive method
Pred.Fk.Naive = function(f,m,k.show=3){
  kmax = length(f)
  est.fk = rep(0,k.show)
  n = sum(f*(1:kmax))
  est.p=1:length(f)
  est.p = est.p/n
  for(j in 1:k.show){
    # est.fk[j] = choose(m,j)*sum(f*est.p^j*(1-est.p)^(n+m-j))
    est.fk[j] = sum(f*dbinom(j,m,est.p)*(1-est.p)^n)
  }
  est.fk
}



##### the unweighted method
Pred.Fk.unweighted = function(f,m,b,f0, k.show=3){
  kmax = length(f)
  
  cut.pts = max(10,kmax)
  
  n = sum(f*(1:kmax))
  d = rep(0,kmax) 
  ## estimated fk
  est.f = rep(0,k.show) 
  
  
  if(f[2]==0){
    f1 <- max(f[1] - 1, 0)
    f2 <- 1
  } else {
    f1 = f[1]
    f2 = f[2]
  }
  d0 = f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*max(f2,1)))/f0
  

  ### relative abundance estimation by Chao et al. (2015)
  ind = 1:kmax
  d = ind/n*(1-b[1]*exp(-b[2]*ind))
  
  for(k in 1:k.show){
  #     est.f[k] = choose(m,k)*sum(
  #       f*d^k*(1-d)^(n+m-k)
  #     );
  #   est.f[k] = est.f[k] + choose(m,k)*f0*d0^k*(1-d0)^(n+m-k)

    est.f[k] = sum(f*dbinom(k,m,d)*(1-d)^n)+f0*dbinom(k,m,d0)*(1-d0)^n;
  }# loop: k
  
  return(est.f)
}


##### main function 
## m: the total number of individuals in an additional sample
## f: species frequency counts data
## xi: species abundance data
Pred.abundance.rare = function(boot.rep = 100,f=NULL,xi=NULL,m,k.show = 3){
  if(is.null(f)*is.null(xi)) {
    print("Please input either species frequency counts data or species abundance data!!")
    return(NULL);
  }
  if(is.null(f)) {
    Xi = xi
    f = X.to.f(Xi)
  }
  if(is.null(xi)) {
    # f=Q
    Xi = f.to.X(f)
  }
  # #### Point estimates  
  n=sum(Xi)
  a.p.1 = DetAbu(Xi, zero=FALSE)
  est.R = Pred.Fk.BW(f=f, m=m, b=a.p.1, k.show)
  est.f0 = SpEst.Chao1.abun(f)-sum(f)
  est.R2 = Pred.Fk.unweighted(f=f, m=m, b=a.p.1, est.f0,k.show)
  est.Naive = Pred.Fk.Naive(f=f,m=m,k.show)

  #### Calculating bootstrap SEs and CIs
  #### Create a space for storing bootstrap samples
  boot.output = NULL
  for(i in 1:boot.rep){
    b.f = boot.abundance.fun(S.hat=est.f0+sum(f), f, b=a.p.1)
    b.Xi = f.to.X(b.f)
    b.a.p.1 = DetAbu(b.Xi, zero=FALSE)
    
    b.n = sum((1:length(b.f))*b.f)
    ## Naive estimator
    b.est.p=1:length(b.f)
    b.est.p = b.est.p/b.n
    b.naive = NULL
    for(j in 1:k.show) {
      b.naive[j] = choose(m,j)*sum(b.f*b.est.p^j*(1-b.est.p)^(b.n+m-j))
    }## loop: j
    
    b.est.f0 = SpEst.Chao1.abun(b.f)-sum(b.f)
    
    b.pred.fk.BW = Pred.Fk.BW(b.f, m, b=b.a.p.1)
    
    b.pred.fk.unweighted = Pred.Fk.unweighted(b.f, m, b=b.a.p.1, b.est.f0)
    
    boot.output = rbind(boot.output,c(proposed=b.pred.fk.BW[1:k.show], naive=b.naive[1:k.show], unweighted=b.pred.fk.unweighted[1:k.show]))
  }### loop: boot.rep
  
  point.est = cbind(proposed=est.R,naive=est.Naive,unweighted=est.R2)
  boot.sd = apply(boot.output, 2, sd, na.rm=T)
  boot.sd = matrix(boot.sd, ncol=3)
  boot.ci = apply(boot.output, 2, quantile, probs=c(0.025,0.975))
  
  Normal.ci = cbind(point.est-1.96*boot.sd,point.est+1.96*boot.sd)
  for(i in 1:nrow(Normal.ci)){
    for(j in 1:ncol(Normal.ci)){
      if(Normal.ci[i,j]<0) Normal.ci[i,j]=0;
    }
  }
  
  output = list()
  output[["Data information"]]=c("Original sample size (n)"=n,"Additional sample size (m)"=m)
  output[["Naive method"]] = round(cbind(1:k.show,point.est[,2],boot.sd[,2],Normal.ci[,c(2,5)]),1)
  output[["Bayesian-weight method"]] = round(cbind(1:k.show,point.est[,1],boot.sd[,1],Normal.ci[,c(1,4)]),1)
  output[["Unweighted method"]] = round(cbind(1:k.show,point.est[,3],boot.sd[,3],Normal.ci[,c(3,6)]),1)

  colnames(output[["Naive method"]])=
    colnames(output[["Bayesian-weight method"]])=
    colnames(output[["Unweighted method"]])=
    c("k","Estimate","Estimated SE","95% lower limit","95% upper limit")
  
   output
} ### end of R function "Pred.incidence.rare"     



################################################
#For incidence-based data
### Chao et al. (2015)
#' DetInc(y, zero=FALSE) is a function of estimating detected species incidence probability.
#' @param y a vector of species incidence frequency
#' @param zero reserves zero frequency or not. Default is FALSE.
#' @return a numerical vector
DetInc <- function(y,nT, zero=FALSE){
  # y <- unlist(y)
  # nT <- max(y)  
  # y <- y[-1]
  Q1 <- sum(y==1)
  Q2 <- sum(y==2)
  Q3 <- sum(y==3)
  if(Q2==0){
    Q1 <- max(Q1 - 1, 0)
    Q2 <- 1
  }
  A1 <- Q1 / nT * ((nT-1)*Q1 / ((nT-1)*Q1 + 2*max(Q2,1)))
  A2 <- Q2 / choose(nT, 2) * ((nT-2)*Q2/((nT-2)*Q2 + 3*max(Q3,1)))^2
  if(zero==FALSE) y <- y[y>0]
  q.solve <- function(q){
    e <- A1 / sum(y/nT*exp(-q*y))
    out <- sum((y/nT * (1 - e * exp(-q*y)))^2) - sum(choose(y,2)/choose(nT,2)) + A2
    abs(out)
  }
  q <- tryCatch(optimize(q.solve, c(0,1))$min, error = function(e) {1})
  e <- A1 / sum(y/nT*exp(-q*y))
  o <- y/nT * (1 - e * exp(-q*y))
  Q0.hat <- ceiling(ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2))  #Chao2 estimator
  return(c(e,q,Q0.hat))
  # o
}


##### Unweighted
Pred.Qk.unweighted = function(Q,nT,u,b,Q0,k.show = 3){
  f0=Q0
  n = nT
  f=Q
  m=u
  kmax = length(f)
  
  cut.pts = max(10,kmax)
  
  ## modified f
  a.f = rep(0,cut.pts)
  ## delta i
  d = rep(0,cut.pts) 
  ## estimated fk
  est.f = rep(0,k.show) 
  
  if(f[2]==0){
    f1 <- max(f[1] - 1, 0)
    f2 <- 1
  } else {
    f1 = f[1]
    f2 = f[2]
  }
  d0 = f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*max(f2,1)))/Q0

  ### relative abundance estimation by Chao et al. (2015)
  ind = 1:kmax
  d = ind/nT*(1-b[1]*exp(-b[2]*ind))
  
  for(k in 1:k.show){
    # est.f[k] = choose(m,k)*sum(f*d^k*(1-d)^(n+m-k));
    # est.f[k] = est.f[k] + choose(m,k)*f0*d0^k*(1-d0)^(n+m-k)

    est.f[k] = sum(f*dbinom(k,m,d)*(1-d)^n)+f0*dbinom(k,m,d0)*(1-d0)^n;
  }# loop: k
  
  return(est.f)
}



##### Proposed Bayesian-weight method
# Pred.Qk.BW = function(Q,nT,u,b){
Pred.Qk.BW = function(Q,nT,u,b,k.show = 3){
  n = nT
  f=Q
  m=u
  kmax = length(f)
  
  cut.pts = max(10,kmax)
  
  ## delta i
  d = rep(0,cut.pts) 
  ## estimated fk
  est.f = rep(0,k.show) 
  
  ### relative abundance estimation by Chao et al. (2015)
  ind = 1:kmax
  d = ind/nT*(1-b[1]*exp(-b[2]*ind))
  
  
  for(i in 1:k.show){## i for k
    for(j in 1:min(i,kmax)){## j for q
      # est.f[i] = est.f[i] +
      #   f[j]*choose(m,i)/(choose(m+n,i)-choose(m,i))*
      #   dbinom(i-j,m,d[j]);
      
      tmp = prod(1+n/(m-(0:(i-1))))
      est.f[i] = est.f[i] +
        f[j]*1/(tmp-1)*dbinom(i-j,m,d[j]);
    }
  }
  
  return(est.f)
}




### Naive method
Pred.Qk.Naive = function(nT,u,f,k.show = 3){
  m.i=u
  n.i=nT
  est.fk = c()
  for(j in 1:k.show) {
    est.p=1:length(f)
    est.p = est.p/nT
    # est.fk[j] = sum(f*(1-est.p)^n.i*choose(m.i,j)*est.p^j*(1-est.p)^(m.i-j))
    
    est.fk[j] = sum(f*dbinom(j,m.i,est.p)*(1-est.p)^n.i)
    
  }
  est.fk
}



### incidence-based Chao et al. (2013): 
#Generating species counts by the bootstrapping method
boot.incidence.fun = function(S.hat, nT, Q, b) {
  f = Q
  D = sum(f)
  kmax = length(f)
  ind = 1:kmax
  d = rep(0,kmax)
  ### the total number of quadrats in the original sample
  n = nT
  
  ### the estimated species richness in the original sample
  S.hat = ceiling(S.hat)
  p0=0;f0=0;
  
  if(f[2]==0){
    f1 <- max(f[1] - 1, 0)
    f2 <- 1
  } else {
    f1 = f[1]
    f2 = f[2]
  }
  #### the estimated sample coverage
  C = 1-f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*max(f2,1)))
  
  d = ind/nT*(1-b[1]*exp(-b[2]*ind))

  cell.p = rep(d, times=f)
  
  if(S.hat>D) {
    f0 = S.hat - D
    p0=(1-C)/f0
    det.p = c(rep(p0,f0),cell.p) 
  } else {
    det.p = cell.p
  }
  sp.count = rep(0,length(det.p))
  
  for(i in 1:length(det.p)){
    sp.count[i] = rbinom(1,n,det.p[i])
  }

  ### bootstrapping sample represented by f1, f2, ....
  f.count = X.to.f(sp.count)
  
  return(f.count)
}



##### main function 
Pred.incidence.rare = function(boot.rep = 100,Q=NULL,xi=NULL,nT,u,k.show = 3){
  ### the number of quadrats in an additional sample
  m = u
  if(is.null(Q)*is.null(xi)) {
    print("Please input either frequency counts data or incidence counts data!!")
    return(NULL);
  }
  if(is.null(Q)) {
    Xi = xi
    f = X.to.f(Xi)
  }
  if(is.null(xi)) {
    f=Q
    Xi = f.to.X(f)
  }
  #### Point estimates  
  out.chao.pi.q = DetInc(Xi,nT, zero=FALSE)
  a.p.1 = out.chao.pi.q[1:2]
  Q0.est = out.chao.pi.q[3]
  
  est.R = Pred.Qk.BW(f,nT,m,b=a.p.1,k.show)
  est.R2 = Pred.Qk.unweighted(f,nT,m,b=a.p.1,Q0.est,k.show)
  est.Naive = Pred.Qk.Naive(nT,m,f,k.show)
  
  #### Calculating bootstrap SEs and CIs
  #### Create a space for storing bootstrap samples
  boot.output = NULL
  
  for(i in 1:boot.rep){
    b.f = boot.incidence.fun(S.hat = sum(f)+Q0.est, nT, f, b=a.p.1)
    b.Xi = f.to.X(b.f)
    
    b.out.chao.pi.q = DetInc(Xi,nT, zero=FALSE)
    b.a.p.1 = b.out.chao.pi.q[1:2]
    b.Q0.est = b.out.chao.pi.q[3]
    
    b.a.p.1 = DetInc(b.Xi, nT, zero=FALSE)
    
    # b.naive = Pred.Qk.Naive(nT, m, b.f, k.show)
    b.est.R = Pred.Qk.BW(b.f,nT,m,b=b.a.p.1, k.show)
    b.est.R2 = Pred.Qk.unweighted(b.f, nT, m, b=b.a.p.1, b.Q0.est, k.show)
    b.est.Naive = Pred.Qk.Naive(nT, m, b.f, k.show)
    
    boot.output = rbind(boot.output,
                        c(proposed=b.est.R[1:k.show], 
                          naive=b.est.Naive[1:k.show], 
                          unweighted=b.est.R2[1:k.show]))
  }### loop: boot.rep
  
  point.est = cbind(proposed=est.R,naive=est.Naive,unweighted=est.R2)
  boot.sd = apply(boot.output, 2, sd, na.rm=T)
  boot.sd = matrix(boot.sd, ncol=3)
  boot.ci = apply(boot.output, 2, quantile, probs=c(0.025,0.975),na.rm=T)
  
  Normal.ci = cbind(point.est-1.96*boot.sd, point.est+1.96*boot.sd)
  # Normal.ci = mapply(max,cbind(point.est-1.96*boot.sd, point.est+1.96*boot.sd),0)
  
  for(i in 1:nrow(Normal.ci)){
    for(j in 1:ncol(Normal.ci)){
      if(Normal.ci[i,j]<0) Normal.ci[i,j]=0; 
    }
  }
  
  output = list()
  output[["Data information"]]=c("Original sample size (t)"=nT,"Additional sample size (u)"=u)
  output[["Naive method"]] = round(cbind(1:k.show,point.est[,2],boot.sd[,2],Normal.ci[,c(2,5)]),1)
  output[["Bayesian-weight method"]] = round(cbind(1:k.show,point.est[,1],boot.sd[,1],Normal.ci[,c(1,4)]),1)
  output[["Unweighted method"]] = round(cbind(1:k.show,point.est[,3],boot.sd[,3],Normal.ci[,c(3,6)]),1)
  
  colnames(output[["Naive method"]])=
    colnames(output[["Bayesian-weight method"]])=
    colnames(output[["Unweighted method"]])=
    c("k","Estimate","Estimated SE","95% lower limit","95% upper limit")
  
  output
} ### end of R function "Pred.incidence.rare"     
#