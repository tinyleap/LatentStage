#' Preselected changepoint, Poisson parameter t.s.
#'
#' This function implements a Viterbi algorithm, given m timestamps, single parameter theta, and n x m emission matrix.
#'
#' @param ts size m vector of timestamps associated with observed outcomes
#' @param theta single parameter, used with ts to inform transition matrix
#' @param obs n x m emission matrix of data observed at different timestamps (column) given different states (row)
#'
#' @return a vector representing the Viterbi path
#'
#' @examples
#' # TO BE DETERMINED
#'
#' @export
preselect.poisson.calc = function(data, pre, minseglen) {
  if (pre >= length(data) | pre <= 0) {
    stop("Changepoint must be located between endpoints of variable range.")
  }
  singledim=function(data, pre, minseglen){
    n=length(data)
    y=c(0,cumsum(data))
    if(y[n+1]==0){
      null=Inf
    }
    else{
      null=2*y[n+1]*log(n) - 2*y[n+1]*log(y[n+1])
    }
    taustar=minseglen:(n-minseglen)
    tmp=2*log(taustar)*y[taustar+1] -2*y[taustar+1]*log(y[taustar+1]) + 2*log(n-taustar)*(y[n+1]-y[taustar+1])-2*(y[n+1]-y[taustar+1])*log((y[n+1]-y[taustar+1]))
    if(sum(is.na(tmp))!=0){
      tmp[which(is.na(tmp))]=Inf
    }
    tau=pre
    taulike=tmp[pre-minseglen+1] # correcting for the fact that we are starting at minseglen
    out=c(tau,null,taulike)
    names(out)=c('cpt','null','alt')
    return(out)
  }
  if(is.null(dim(data))==TRUE){
    cpt=singledim(data,pre,minseglen)
    return(cpt)
  }
  else{
    rep=nrow(data)
    n=ncol(data)
    cpt=matrix(0,ncol=3,nrow=rep)
    for(i in 1:rep){
      cpt[i,]=singledim(data[i,],pre,minseglen)
    }
    colnames(cpt)=c('cpt','null','alt')
    return(cpt)
  }
}

#' @export
preselect.norm.mean.calc = function(data,pre,minseglen){
  if (pre >= length(data) | pre <= 0) {
    stop("Changepoint must be located between endpoints of variable range.")
  }
  singledim=function(data,pre,minseglen){
    n=length(data)
    y=c(0,cumsum(data))
    y2=c(0,cumsum(data^2))
    null=y2[n+1]-y[n+1]^2/n
    taustar=minseglen:(n-minseglen+1)
    tmp=y2[taustar+1]-y[taustar+1]^2/taustar + (y2[n+1]-y2[taustar+1]) - ((y[n+1]-y[taustar+1])^2)/(n-taustar)
    tau=pre
    taulike=tmp[pre-minseglen+1] # correcting for the fact that we are starting at minseglen
    out=c(tau,null,taulike)
    names(out)=c('cpt','null','alt')
    return(out)
  }
  if(is.null(dim(data))==TRUE){
    cpt=singledim(data,pre,minseglen)
    return(cpt)
  } else{
    rep=nrow(data)
    n=ncol(data)
    cpt=matrix(0,ncol=3,nrow=rep)
    for(i in 1:rep){
      cpt[i,]=singledim(data[i,],pre,minseglen)
    }
    colnames(cpt)=c('cpt','null','alt')
    return(cpt)
  }
}

#' @export
preselect.norm.var.calc = function(data,pre,know.mean,mu,minseglen){
  if (pre >= length(data) | pre <= 0) {
    stop("Changepoint must be located between endpoints of variable range.")
  }
  if(is.null(dim(data))==TRUE){
    if((know.mean==FALSE)&(is.na(mu))){
      mu=mean(coredata(data))
    }
  } else {
    rep=nrow(data)
    if(length(mu)!=rep){
      mu=rep(mu,rep)
    }
    for(i in 1:rep){
      if((know.mean==FALSE)&(is.na(mu[i]))){
        mu=mean(coredata(data[i,]))
      }
    }
  }
  n=length(data)
  y=c(0,cumsum((data-mu)^2))
  null=n*log(y[n+1]/n)
  taustar=minseglen:(n-minseglen+1)
  sigma1=y[taustar+1]/taustar
  neg=sigma1<=0
  sigma1[neg==TRUE]=1*10^(-10)
  sigman=(y[n+1]-y[taustar+1])/(n-taustar)
  neg=sigman<=0
  sigman[neg==TRUE]=1*10^(-10)
  tmp=taustar*log(sigma1) + (n-taustar)*log(sigman)
  tau=pre
  taulike=tmp[pre-minseglen+1]
  out=c(tau,null,taulike)
  names(out)=c('cpt','null','alt')
  return(out)
}

#' @export
preselect.norm.meanvar.calc = function(data,pre,minseglen){
  singledim=function(data,pre,minseglen){
    n=length(data)
    y=c(0,cumsum(data))
    y2=c(0,cumsum((data)^2))
    null=n*log((y2[n+1]-(y[n+1]^2/n))/n)
    taustar=minseglen:(n-minseglen+1)
    sigma1=((y2[taustar+1]-(y[taustar+1]^2/taustar))/taustar)
    neg=sigma1<=0
    sigma1[neg==TRUE]=1*10^(-10)
    sigman=((y2[n+1]-y2[taustar+1])-((y[n+1]-y[taustar+1])^2/(n-taustar)))/(n-taustar)
    neg=sigman<=0
    sigman[neg==TRUE]=1*10^(-10)
    tmp=taustar*log(sigma1) + (n-taustar)*log(sigman)
    tau=pre
    taulike=tmp[pre-minseglen+1] # correcting for the fact that we are starting at minseglen
    out=c(tau,null,taulike)
    names(out)=c('cpt','null','alt')
    return(out)
  }
  if(is.null(dim(data))==TRUE){
    cpt=singledim(data,pre,minseglen)
    return(cpt)
  }
  else{
    rep=nrow(data)
    n=ncol(data)
    cpt=matrix(0,ncol=3,nrow=rep)
    for(i in 1:rep){
      cpt[i,]=singledim(data[i,],pre,minseglen)
    }
    colnames(cpt)=c('cpt','null','alt')
    return(cpt)
  }
}

#' p-Value calculation for binary model
#' 
#' DESC - Reminder to myself: x and y need to be sorted beforehand
#' 
#' DETAILS - Also, create function or text for separating data to find multiple
#' changepoints.
#' 
#' @export
p.value <- function(y, x, cpt) {
  N = length(y)
  n = sum(y)
  num = (n/N)^n * ((N-n)/N)^(N-n)
  chex = min(which(x==cpt))
  N_1 = chex-1
  N_2 = N - chex
  n_1 = sum(y[1:(chex-1)])
  n_2 = sum(y[chex:N])
  den = (n_1/N_1)^n_1 * ((N_1-n_1)/N_1)^(N_1-n_1) * (n_2/N_2)^n_2 * ((N_2-n_2)/N_2)^(N_2-n_2)
  lambda = num / den
  ts = -2*log(lambda)
  p = 1 - pchisq(ts, df=2)
  return(p)
}

#' @export
change.normal = function(v, pre=NA, type="MeanVar", know.mean=FALSE, mu=NA, minseglen=2) {
  if(is.null(dim(v))==TRUE){
    n=length(v)
  } else{
    n=ncol(v)
  }
  if(n<4){stop('Data must have at least 4 observations to fit a changepoint model.')}
  if(n<(2*minseglen)){stop('Minimum segment length is too large to include a change in this data')}
  if (is.na(pre) == FALSE) {
    if (type=="Mean") {
      cpt.pen = penalty_decision(penalty="MBIC", pen.value=0, n=n, diffparam=1, asymcheck="mean.norm", method="AMOC")
      if(is.null(dim(v))==TRUE){ # single dataset
        tmp=preselect.norm.mean.calc(coredata(v),pre,minseglen)
        tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
        alogn=(2*log(log(n)))^(-(1/2))
        blogn=(alogn^(-1))+(1/2)*alogn*log(log(log(n)))  # Chen & Gupta (2000) pg10
        p.value=1-exp(-2*(pi^(1/2))*exp(-alogn*sqrt(abs(tmp[2]-tmp[3]))+(alogn^{-1})*blogn))-exp(-2*(pi^(1/2))*exp((alogn^{-1})*blogn))
        my_results = c("CP Location", "CP Penalty (AMOC)", "CP p-Value", "MBIC Penalty Threshold")
        my_row = c(tmp[1], tmp[2] - tmp[3], p.value, cpt.pen)
        my_results = rbind(my_results, my_row)
        return(my_results)
      } else{ 
        tmp=preselect.norm.mean.calc(v,pre,minseglen)
        tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
        alogn=(2*log(log(n)))^(-(1/2))
        blogn=(alogn^(-1))+(1/2)*alogn*log(log(log(n)))  # Chen & Gupta (2000) pg10
        p.value = 1-exp(-2*(pi^(1/2))*exp(-alogn*sqrt(abs(tmp[,2]-tmp[,3]))+(alogn^{-1})*blogn))-exp(-2*(pi^(1/2))*exp((alogn^{-1})*blogn))
        rep=nrow(v)
        my_results = c("Changepoint #", "CP Location", "CP Penalty (AMOC)", "CP p-Value", "MBIC Penalty Threshold")
        my_row=list()
        for(i in 1:rep){
          my_row[[i]] = c(i, tmp[i,1], tmp[i,2] - tmp[i,3], p.value[i], cpt.pen)
        }
        my_results = rbind(my_results, my_row)
        return(my_results)
      }
    } else if (type=="Var") {
      mu=mu[1]
      cpt.pen = penalty_decision(penalty="MBIC", pen.value=0, n=n, diffparam=1, asymcheck="var.norm", method="AMOC")
      if(is.null(dim(v))==TRUE){
        if((know.mean==FALSE)&(is.na(mu))){
          mu=mean(coredata(v))
        }
        tmp=preselect.norm.var.calc(coredata(v),pre,know.mean,mu,minseglen)
        tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
        alogn=sqrt(2*log(log(n)))
        blogn=2*log(log(n))+ (log(log(log(n))))/2 - log(gamma(1/2))
        p.value=1-exp(-2*exp(-alogn*sqrt(abs(tmp[2]-tmp[3]))+blogn))-exp(-2*exp(blogn))  # Chen & Gupta (2000) pg27
        my_results = c("CP Location", "CP Penalty (AMOC)", "CP p-Value", "MBIC Penalty Threshold")
        my_row = c(tmp[1], tmp[2] - tmp[3], p.value, cpt.pen)
        my_results = rbind(my_results, my_row)
        return(my_results)
      } else{ 
        rep=nrow(v)
        tmp=matrix(0,ncol=3,nrow=rep)
        if(length(mu)!=rep){
          mu=rep(mu,rep)
        }
        for(i in 1:rep){
          if((know.mean==FALSE)&(is.na(mu[i]))){
            mu=mean(coredata(v[i,]))
          }
          tmp[i,]=preselect.norm.var.calc(v[i,],pre,know.mean,mu[i],minseglen)
        }
        tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
        alogn=sqrt(2*log(log(n)))
        blogn=2*log(log(n))+ (log(log(log(n))))/2 - log(gamma(1/2))
        p.value=1-exp(-2*exp(-alogn*sqrt(abs(tmp[,2]-tmp[,3]))+blogn))-exp(-2*exp(blogn))  # Chen & Gupta (2000) pg27
        my_results = c("Changepoint #", "CP Location", "CP Penalty (AMOC)", "CP p-Value", "MBIC Penalty Threshold")
        my_row=list()
        for(i in 1:rep){
          my_row[[i]] = c(i, tmp[i,1], tmp[i,2] - tmp[i,3], p.value[i], cpt.pen)
        }
        my_results = rbind(my_results, my_row)
        return(my_results)
      }
    } else if (type=="MeanVar") {
      cpt.pen = penalty_decision(penalty="MBIC", pen.value=0, n=n, diffparam=2, asymcheck="meanvar.norm", method="AMOC")
      if(is.null(dim(v))==TRUE){
        tmp=preselect.norm.meanvar.calc(coredata(v),pre,minseglen)
        tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
        alogn=sqrt(2*log(log(n)))
        blogn=2*log(log(n))+ (log(log(log(n))))
        p.value=1-exp(-2*exp(-alogn*sqrt(abs(tmp[2]-tmp[3]))+blogn))-exp(-2*exp(blogn))   # Chen & Gupta (2000) pg54
        my_results = c("CP Location", "CP Penalty (AMOC)", "CP p-Value", "MBIC Penalty Threshold")
        my_row = c(tmp[1], tmp[2] - tmp[3], p.value, cpt.pen)
        my_results = rbind(my_results, my_row)
        return(my_results)
      } else{ 
        tmp=preselect.norm.meanvar.calc(v,pre,minseglen)
        tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
        alogn=sqrt(2*log(log(n)))
        blogn=2*log(log(n))+ (log(log(log(n))))
        p.value=1-exp(-2*exp(-alogn*sqrt(abs(tmp[,2]-tmp[,3]))+blogn))-exp(-2*exp(blogn))   # Chen & Gupta (2000) pg54
        rep=nrow(v)
        my_results = c("Changepoint #", "CP Location", "CP Penalty (AMOC)", "CP p-Value", "MBIC Penalty Threshold")
        my_row=list()
        for(i in 1:rep){
          my_row[[i]] = c(i, tmp[i,1], tmp[i,2] - tmp[i,3], p.value[i], cpt.pen)
        }
        my_results = rbind(my_results, my_row)
        return(my_results)
      }
    } else {
      stop("Incorrect change type. Please select from 'Mean', 'Var', or 'MeanVar'.")
    }
  } else {
    if (type=="Mean") {
      cpt = cpt.mean(v, test.stat="Normal", method="BinSeg", Q = n / 2 + 1, minseglen=minseglen)
    } else if (type=="Var") {
      cpt = cpt.var(v, test.stat="Normal", method="BinSeg", know.mean = know.mean, mu=mu, Q = n / 2 + 1, minseglen=minseglen)
    } else if (type=="MeanVar") {

      cpt = cpt.meanvar(v, test.stat="Normal", method="BinSeg", Q = n / 2 + 1, minseglen=minseglen)
    } else {
      stop("Incorrect change type. Please select from 'Mean', 'Var', or 'MeanVar'.")
    }
    my_results = c("Changepoint #", "CP Location", "CP Penalty (BinSeg)", "MBIC Penalty Threshold")
    for (i in 1:(n/2 + 1)) {
      if (cpt@pen.value.full[i] >= cpt@pen.value) {
        my_row = c(i, cpt@cpts.full[i,i], cpt@pen.value.full[i], cpt@pen.value)
        my_results = rbind(my_results, my_row)
      }
    }
    return(my_results)
  }
}

#' @export
change.binary = function(formula.1, formula.2, df) {
  cpt1 = chngptm(formula.1 = formula.1, formula.2=formula.2, df, family="binomial", type="step")
  p.value = p.value(formula.1, formula.2, cpt1$coefficients)
  my_results = c(names(cpt1$coefficients), "CP p-Value")
  my_row = c(cpt1$coefficients, cpt$p.value)
  my_results = rbind(my_results, my_row)
  return(my_results)
}

#' @export
change.poisson = function(v, pre=NA, minseglen=2) {
  if((sum(v<0)>0)){stop('Poisson test statistic requires positive data')}
  if(sum(as.integer(v)==v)!=length(v)){stop('Poisson test statistic requires integer data')}
  if(is.null(dim(v))==TRUE){
    n=length(v)
  }
  else{
    n=ncol(v)
  }
  if(n<4){stop('Data must have at least 4 observations to fit a changepoint model.')}
  if(n<(2*minseglen)){stop('Minimum segment length is too large to include a change in this data')}
  if (is.na(pre) == FALSE) {
    cpt.pen = penalty_decision(penalty="MBIC", pen.value=0, n=n, diffparam=1, asymcheck="meanvar.poisson", method="AMOC")
    if(is.null(dim(v))==TRUE){
      tmp=preselect.poisson.calc(coredata(v),pre,minseglen)
      tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
      my_results = c("CP Location", "CP Penalty (AMOC)", "MBIC Penalty Threshold")
      my_row = c(tmp[1], tmp[2] - tmp[3], cpt.pen)
      my_results = rbind(my_results, my_row)
      return(my_results)
    }
    else{ 
      tmp=preselect.poisson.calc(v,pre,minseglen)
      tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
      rep=nrow(v)
      my_results = c("Changepoint #", "CP Location", "CP Penalty (AMOC)", "MBIC Penalty Threshold")
      my_row=list()
      for(i in 1:rep){
        my_row[[i]]== c(i, tmp[i,1], tmp[i,2] - tmp[i,3], cpt.pen)
      }
      my_results = rbind(my_results, my_row)
      return(my_results)
    }
  } else {
    cpt = cpt.meanvar(v, test.stat="Poisson", method="BinSeg", minseglen=minseglen, Q= n / 2 + 1)
    my_results = c("Changepoint #", "CP Location", "CP Penalty (BinSeg)", "MBIC Penalty Threshold")
    for (i in 1:(n/2 + 1)) {
      if (cpt@pen.value.full[i] >= cpt@pen.value) {
        my_row = c(i, cpt@cpts.full[i,i], cpt@pen.value.full[i], cpt@pen.value)
        my_results = rbind(my_results, my_row)
      }
    }
    return(my_results)
  }
}