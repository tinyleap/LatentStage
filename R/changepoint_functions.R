#' Preselect changepoint for Poisson rate parameter
#'
#' Used in \code{\link{change.poisson}} to calculate significance of a
#' user-selected location as a changepoint.
#'
#' The code behind calculating the significance of a changepoint is derived from
#' \code{\link{changepoint}} package. Insert reference to original paper here,
#' specifically on BIC for changepoint models.
#'
#' @param data The data on which you wish to test the significance of a
#'   preselected changepoint.
#' @param pre The changepoint in the data for which you wish to test the
#'   significance.
#'
#' @return A data frame of values needed to calculate significance of location
#'   as a changepoint.
#'
#' @export
preselect.poisson.calc = function(data, pre) {
  if (pre >= length(data) | pre <= 0) {
    stop("Changepoint must be located between endpoints of variable range.")
  }
  singledim=function(data, pre){
    n=length(data)
    y=c(0,cumsum(data))
    if(y[n+1]==0){
      null=Inf
    }
    else{
      null=2*y[n+1]*log(n) - 2*y[n+1]*log(y[n+1])
    }
    taustar=2:(n-2)
    tmp=2*log(taustar)*y[taustar+1] -2*y[taustar+1]*log(y[taustar+1]) + 2*log(n-taustar)*(y[n+1]-y[taustar+1])-2*(y[n+1]-y[taustar+1])*log((y[n+1]-y[taustar+1]))
    if(sum(is.na(tmp))!=0){
      tmp[which(is.na(tmp))]=Inf
    }
    tau=pre
    taulike=tmp[pre-2+1]
    out=c(tau,null,taulike)
    names(out)=c('cpt','null','alt')
    return(out)
  }
  if(is.null(dim(data))==TRUE){
    cpt=singledim(data,pre)
    return(cpt)
  }
  else{
    rep=nrow(data)
    n=ncol(data)
    cpt=matrix(0,ncol=3,nrow=rep)
    for(i in 1:rep){
      cpt[i,]=singledim(data[i,],pre)
    }
    colnames(cpt)=c('cpt','null','alt')
    return(cpt)
  }
}

#' Preselect changepoint for Normal mean
#'
#' Used in \code{\link{change.normal}} to calculate significance of a
#' user-selected location as a changepoint. Specific to \code{type = "Mean"}.
#'
#' The code behind calculating the significance of a changepoint is derived from
#' \code{\link{changepoint}} package. Insert reference to original paper here,
#' specifically on BIC for changepoint models.
#'
#' @param data The data on which you wish to test the significance of a
#'   preselected changepoint.
#' @param pre The changepoint in the data for which you wish to test the
#'   significance.
#'
#' @return A data frame of values needed to calculate significance of location
#'   as a changepoint.
#'
#' @export
preselect.norm.mean.calc = function(data,pre){
  if (pre >= length(data) | pre <= 0) {
    stop("Changepoint must be located between endpoints of variable range.")
  }
  singledim=function(data,pre){
    n=length(data)
    y=c(0,cumsum(data))
    y2=c(0,cumsum(data^2))
    null=y2[n+1]-y[n+1]^2/n
    taustar=2:(n-2+1)
    tmp=y2[taustar+1]-y[taustar+1]^2/taustar + (y2[n+1]-y2[taustar+1]) - ((y[n+1]-y[taustar+1])^2)/(n-taustar)
    tau=pre
    taulike=tmp[pre-2+1]
    out=c(tau,null,taulike)
    names(out)=c('cpt','null','alt')
    return(out)
  }
  if(is.null(dim(data))==TRUE){
    cpt=singledim(data,pre)
    return(cpt)
  } else{
    rep=nrow(data)
    n=ncol(data)
    cpt=matrix(0,ncol=3,nrow=rep)
    for(i in 1:rep){
      cpt[i,]=singledim(data[i,],pre)
    }
    colnames(cpt)=c('cpt','null','alt')
    return(cpt)
  }
}

#' Preselect changepoint for Normal variance
#'
#' Used in \code{\link{change.normal}} to calculate significance of a
#' user-selected location as a changepoint. Specific to \code{type = "Var"}.
#'
#' The code behind calculating the significance of a changepoint is derived from
#' \code{\link{changepoint}} package. Insert reference to original paper here,
#' specifically on BIC for changepoint models.
#'
#' @param data The data on which you wish to test the significance of a
#'   preselected changepoint.
#' @param pre The changepoint in the data for which you wish to test the
#'   significance.
#' @param know.mean Taken from \code{\link{change.normal}}, Boolean variable
#'   indicating whether or not you wish to set value for "true mean" of data.
#' @param mu Taken from \code{\link{change.normal}}, if \code{know.mean = TRUE},
#'   single value that represents "true mean" of data.
#'
#' @return A data frame of values needed to calculate significance of location
#'   as a changepoint.
#'
#' @export
preselect.norm.var.calc = function(data,pre,know.mean,mu){
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
  taustar=2:(n-2+1)
  sigma1=y[taustar+1]/taustar
  neg=sigma1<=0
  sigma1[neg==TRUE]=1*10^(-10)
  sigman=(y[n+1]-y[taustar+1])/(n-taustar)
  neg=sigman<=0
  sigman[neg==TRUE]=1*10^(-10)
  tmp=taustar*log(sigma1) + (n-taustar)*log(sigman)
  tau=pre
  taulike=tmp[pre-2+1]
  out=c(tau,null,taulike)
  names(out)=c('cpt','null','alt')
  return(out)
}

#' Preselect changepoint for Normal mean and variance
#'
#' Used in \code{\link{change.normal}} to calculate significance of a
#' user-selected location as a changepoint. Specific to \code{type = "MeanVar"}.
#'
#' The code behind calculating the significance of a changepoint is derived from
#' \code{\link{changepoint}} package. Insert reference to original paper here,
#' specifically on BIC for changepoint models.
#'
#' @param data The data on which you wish to test the significance of a
#'   preselected changepoint.
#' @param pre The changepoint in the data for which you wish to test the
#'   significance.
#'
#' @return A data frame of values needed to calculate significance of location
#'   as a changepoint.
#'
#' @export
preselect.norm.meanvar.calc = function(data,pre){
  singledim=function(data,pre){
    n=length(data)
    y=c(0,cumsum(data))
    y2=c(0,cumsum((data)^2))
    null=n*log((y2[n+1]-(y[n+1]^2/n))/n)
    taustar=2:(n-2+1)
    sigma1=((y2[taustar+1]-(y[taustar+1]^2/taustar))/taustar)
    neg=sigma1<=0
    sigma1[neg==TRUE]=1*10^(-10)
    sigman=((y2[n+1]-y2[taustar+1])-((y[n+1]-y[taustar+1])^2/(n-taustar)))/(n-taustar)
    neg=sigman<=0
    sigman[neg==TRUE]=1*10^(-10)
    tmp=taustar*log(sigma1) + (n-taustar)*log(sigman)
    tau=pre
    taulike=tmp[pre-2+1]
    out=c(tau,null,taulike)
    names(out)=c('cpt','null','alt')
    return(out)
  }
  if(is.null(dim(data))==TRUE){
    cpt=singledim(data,pre)
    return(cpt)
  }
  else{
    rep=nrow(data)
    n=ncol(data)
    cpt=matrix(0,ncol=3,nrow=rep)
    for(i in 1:rep){
      cpt[i,]=singledim(data[i,],pre)
    }
    colnames(cpt)=c('cpt','null','alt')
    return(cpt)
  }
}

#' p-Value Calculator for Binary Discrete Data
#'
#' Calculates significance value, i.e. p-value, of a location in data as a
#' changepoint, over range of specified threshold variable.
#'
#' Include information on likelihood ratio test for binary discrete data, as
#' well as plugging the parameter of said test into the chi-squared distribution
#' with degrees of freedom equal to 3.
#'
#' @param v Vector that contains response variable of interest.
#' @param t Vector of same length containing threshold variable.
#' @param cpt Value representing the location in data (either preselected or
#'   determined rigorously) on which to test a changepoint.
#'
#' @return Single value, the p-value of a location as changepoint of y, over t.
#'
#' @export
p.value <- function(v, t, cpt) {
  N = length(v)
  n = sum(v)
  num = (n/N)^n * ((N-n)/N)^(N-n)
  chex = min(which(t==cpt))
  N_1 = chex-1
  N_2 = N - chex
  n_1 = sum(v[1:N_1]) #NOTE: This is not correct! Fix this as soon as possible!
  n_2 = sum(v[(N-N_2):N])
  den = (n_1/N_1)^n_1 * ((N_1-n_1)/N_1)^(N_1-n_1) * (n_2/N_2)^n_2 * ((N_2-n_2)/N_2)^(N_2-n_2)
  lambda = num / den
  ts = -2*log(lambda)
  p = 1 - pchisq(ts, df=2)
  return(p)
}

#' Changepoint Model with Normal Distribution
#'
#' Fit changepoint model on data to test if there is a significant change in
#' mean and/or variance of data.
#'
#' Proceed to describe details of changepoint model, changepoint package,
#' preselect option, what to observe in data that can inform which changes in
#' parameter may be the best fit, and anything else.
#'
#' @param v Vector that contains response variable of interest.
#' @param pre A value  represents preselected changepoint for any model, should
#'   the user wish to test a specific point in data; NA by default.
#' @param type Choose between "Mean", "Var", and "MeanVar", regarding the
#'   parameter(s) for which you wish to test significant change in value;
#'   "MeanVar" by default.
#' @param know.mean If \code{type="Var"}, Boolean variable indicating whether or
#'   not you wish to set value for "true mean" of data; FALSE by default.
#' @param mu If \code{know.mean = TRUE}, single value that represents "true
#'   mean" of data; NA by default.
#'
#' @return List of the following: Data frame that contains location of
#'   changepoint(s) (either preselected or determined rigorously) ranked by
#'   significance and its/their associated significance value(s); Single value
#'   representing the MBIC threshold by which to base significance of a
#'   changepoint (i.e., penalty values LARGER than MBIC threshold indicate
#'   significance)
#'
#' @examples
#' a = rnorm(100,0,25/100)
#' b = rep(0,100)
#' for (i in 1:100) {
#'   b[i] = rnorm(1,0,(i/2+25)/100)
#' }
#' c = rnorm(100,0,75/100)
#' d = rep(0,100)
#' for (i in 1:100)
#'   d[i] = rnorm(1,0,(75-i/2)/100)
#' e = rnorm(100,0,25/100)
#' v = c(a,b,c,d,e)
#'
#' change.normal(v, type="Var", know.mean=TRUE, mu=0)
#' change.normal(v, type="Var", know.mean=TRUE, mu=0, pre=200)
#'
#' @export
change.normal = function(v, pre=NA, type="MeanVar", know.mean=FALSE, mu=NA) {
  if(is.null(dim(v))==TRUE) {
    n=length(v)
  } else{
    n=ncol(v)
  }
  if(n<4){stop('Data must have at least 4 observations to fit a changepoint model.')}
  if(n<(2*2)){stop('Minimum segment length is too large to include a change in this data')}
  if (is.na(pre) == FALSE) {
    if (type=="Mean") {
      cpt.pen = penalty_decision(penalty="MBIC", pen.value=0, n=n, diffparam=1, asymcheck="mean.norm", method="AMOC")
      if(is.null(dim(v))==TRUE) {
        tmp=preselect.norm.mean.calc(coredata(v),pre)
        tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
        alogn=(2*log(log(n)))^(-(1/2))
        blogn=(alogn^(-1))+(1/2)*alogn*log(log(log(n)))  # Chen & Gupta (2000) pg10
        p.value=1-exp(-2*(pi^(1/2))*exp(-alogn*sqrt(abs(tmp[2]-tmp[3]))+(alogn^{-1})*blogn))-exp(-2*(pi^(1/2))*exp((alogn^{-1})*blogn))
        my_row = list(c(tmp[1], tmp[2] - tmp[3], p.value))
        my_df = as.data.frame(do.call(rbind,my_row))
        colnames(my_df) = c("Changepoint Location", "Changepoint Penalty Value")
        return(list(changepoint_data_frame = my_df, MBIC_penalty_threshold = cpt.pen))
      } else{
        tmp=preselect.norm.mean.calc(v,pre)
        tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
        alogn=(2*log(log(n)))^(-(1/2))
        blogn=(alogn^(-1))+(1/2)*alogn*log(log(log(n)))  # Chen & Gupta (2000) pg10
        p.value = 1-exp(-2*(pi^(1/2))*exp(-alogn*sqrt(abs(tmp[,2]-tmp[,3]))+(alogn^{-1})*blogn))-exp(-2*(pi^(1/2))*exp((alogn^{-1})*blogn))
        rep=nrow(v)
        my_row=list()
        for(i in 1:rep){
          my_row[[i]] = c(i, tmp[i,1], tmp[i,2] - tmp[i,3], p.value[i])
        }
        my_df = as.data.frame(do.call(rbind,my_row))
        colnames(my_df) = c("Changepoint #", "Changepoint Location", "Changepoint Penalty Value")
        return(list(changepoint_data_frame = my_df, MBIC_penalty_threshold = cpt.pen))
      }
    } else if (type=="Var") {
      mu=mu[1]
      cpt.pen = penalty_decision(penalty="MBIC", pen.value=0, n=n, diffparam=1, asymcheck="var.norm", method="AMOC")
      if(is.null(dim(v))==TRUE){
        if((know.mean==FALSE)&(is.na(mu))){
          mu=mean(coredata(v))
        }
        tmp=preselect.norm.var.calc(coredata(v),pre,know.mean,mu)
        tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
        alogn=sqrt(2*log(log(n)))
        blogn=2*log(log(n))+ (log(log(log(n))))/2 - log(gamma(1/2))
        p.value=1-exp(-2*exp(-alogn*sqrt(abs(tmp[2]-tmp[3]))+blogn))-exp(-2*exp(blogn))  # Chen & Gupta (2000) pg27
        my_row = list(c(tmp[1], tmp[2] - tmp[3], p.value))
        my_df = as.data.frame(do.call(rbind,my_row))
        colnames(my_df) = c("Changepoint Location", "Changepoint Penalty Value")
        return(list(changepoint_data_frame = my_df, MBIC_penalty_threshold = cpt.pen))
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
          tmp[i,]=preselect.norm.var.calc(v[i,],pre,know.mean,mu[i])
        }
        tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
        alogn=sqrt(2*log(log(n)))
        blogn=2*log(log(n))+ (log(log(log(n))))/2 - log(gamma(1/2))
        p.value=1-exp(-2*exp(-alogn*sqrt(abs(tmp[,2]-tmp[,3]))+blogn))-exp(-2*exp(blogn))
        my_row=list()
        for(i in 1:rep){
          my_row[[i]] = c(i, tmp[i,1], tmp[i,2] - tmp[i,3], p.value[i])
        }
        my_df = as.data.frame(do.call(rbind,my_row))
        colnames(my_df) = c("Changepoint #", "Changepoint Location", "Changepoint Penalty Value")
        return(list(changepoint_data_frame = my_df, MBIC_penalty_threshold = cpt.pen))
      }
    } else if (type=="MeanVar") {
      cpt.pen = penalty_decision(penalty="MBIC", pen.value=0, n=n, diffparam=2, asymcheck="meanvar.norm", method="AMOC")
      if(is.null(dim(v))==TRUE){
        tmp=preselect.norm.meanvar.calc(coredata(v),pre)
        tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
        alogn=sqrt(2*log(log(n)))
        blogn=2*log(log(n))+ (log(log(log(n))))
        p.value=1-exp(-2*exp(-alogn*sqrt(abs(tmp[2]-tmp[3]))+blogn))-exp(-2*exp(blogn))
        my_row = list(c(tmp[1], tmp[2] - tmp[3], p.value))
        my_df = as.data.frame(do.call(rbind,my_row))
        colnames(my_df) = c("Changepoint Location", "Changepoint Penalty Value")
        return(list(changepoint_data_frame = my_df, MBIC_penalty_threshold = cpt.pen))
      } else{
        tmp=preselect.norm.meanvar.calc(v,pre)
        tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
        alogn=sqrt(2*log(log(n)))
        blogn=2*log(log(n))+ (log(log(log(n))))
        p.value=1-exp(-2*exp(-alogn*sqrt(abs(tmp[,2]-tmp[,3]))+blogn))-exp(-2*exp(blogn))
        rep=nrow(v)
        my_row=list()
        for(i in 1:rep){
          my_row[[i]] = c(i, tmp[i,1], tmp[i,2] - tmp[i,3], p.value[i])
        }
        my_df = as.data.frame(do.call(rbind,my_row))
        colnames(my_df) = c("Changepoint #", "Changepoint Location", "Changepoint Penalty Value")
        return(list(changepoint_data_frame = my_df, MBIC_penalty_threshold = cpt.pen))
      }
    } else {
      stop("Incorrect change type. Please select from 'Mean', 'Var', or 'MeanVar'.")
    }
  } else {
    if (type=="Mean") {
      cpt = cpt.mean(v, test.stat="Normal", method="BinSeg", Q = n / 2 + 1, minseglen=2)
    } else if (type=="Var") {
      cpt = cpt.var(v, test.stat="Normal", method="BinSeg", know.mean = know.mean, mu=mu, Q = n / 2 + 1, minseglen=2)
    } else if (type=="MeanVar") {
      cpt = cpt.meanvar(v, test.stat="Normal", method="BinSeg", Q = n / 2 + 1, minseglen=2)
    } else {
      stop("Incorrect change type. Please select from 'Mean', 'Var', or 'MeanVar'.")
    }
    my_df = list()
    j = 1
    for (i in 1:(n/2 + 1)) {
      if (cpt@pen.value.full[i] >= cpt@pen.value) {
        my_row = c(i, cpt@cpts.full[i,i], cpt@pen.value.full[i])
        my_df[[j]] = my_row
        j = j + 1
      }
    }
    if (length(my_df) != 0) {
      my_df = as.data.frame(do.call(rbind,my_df))
      colnames(my_df) = c("Changepoint #", "Changepoint Location", "Changepoint Penalty Value")
      return(list(changepoint_data_frame = my_df, MBIC_penalty_threshold = cpt@pen.value))
    } else {
      return("No significant changepoints were found. If you wish to see the significance of a particular changepoint, please use the 'pre' parameter.")
    }
  }
}

#' Changepoint Model with Binary Discrete Distribution
#'
#' Fit changepoint model on data to test if there is a significant change in
#' binary probability of one outcome versus other outcome.
#'
#' Proceed to describe details of changepoint model, chngpt package,
#' \code{\link{p.value}}, and anything else from Rmd that may matter.
#'
#' @param v Vector that contains response variable of interest.
#' @param pre A value represents preselected changepoint for any model, should
#'   the user wish to test a specific point in data; NA by default.
#' @param t Vector that contains threshold variable over which to study change
#'   (ex. Time, outcome count in succession); default is sequence from 1 to
#'   total number of observations.
#' @param x Vector that contains hidden parameter, should you wish to test an
#'   interaction of original data with another variable (i.e., a hierarchical
#'   model); NA by default.
#'
#' @return Data frame that contains location of changepoint (either preselected
#'   or determined rigorously) and its associated significance value.
#'
#' @examples
#' r = rnorm(400,0,1) * 10
#' r1 = r[1:200]
#' r2 = r[201:400]
#' a = rbinom(200,1,
#'            exp(r1)/(1+exp(r1)))
#' b = rbinom(200,1,
#'            exp(r2+1)/(1+exp(r2+2)))
#' d = data.frame(x=r,y=c(a,b),t=1:400)
#'
#' change.binary(v = d$y, t= d$t, x = d$x)
#' change.binary(v = d$y, t= d$t, pre = 200)
#'
#' @export
change.binary = function(v, t = seq(1, length(v)), x = rep(NA, length(v)), pre=NA) {
  data = data.frame(v = v, t = t, x = x)
  if (is.na(pre)==TRUE) {
    if (all(is.na(x)) == TRUE) {
      cpt1 = chngptm(formula.1 = v ~ 1, formula.2=~t, data, family="binomial", type="step")
    } else {
      cpt1 = chngptm(formula.1 = v ~ x, formula.2=~t, data, family="binomial", type="step")
    }
    p.value = p.value(v, t, cpt1$coefficients[["chngpt"]])
    my_row = list(c(cpt1$coefficients[["chngpt"]], p.value))
    my_df = as.data.frame(do.call(rbind,my_row))
    colnames(my_df) = c("Changepoint Location", "Changepoint p-Value")
    return(my_df)
  } else {
    if (any(is.na(x)) == FALSE) {
      stop("Unfortunately, the preselected changepoint method does not work with a hidden variable x for the Binary statistic.")
    }
    p.value = p.value(v, t, pre)
    my_row = list(c(pre, p.value))
    my_df = as.data.frame(do.call(rbind,my_row))
    colnames(my_df) = c("Changepoint Location", "Changepoint p-Value")
    return(my_df)
  }
}

#' Changepoint Model with Poisson Distribution
#'
#' Fit changepoint model on data to test if there is a significant change in
#' Poisson rate parameter of data.
#'
#' Proceed to describe details of changepoint model, changepoint package,
#' preselect option, what to observe in data that indicates a Poisson
#' changepoint model may be a good fit, and anything else.
#'
#' @param v Vector that contains response variable of interest.
#' @param pre A value  represents preselected changepoint for any model, should
#'   the user wish to test a specific point in data, NA by default.
#'
#' @return List of the following: Data frame that contains location of
#'   changepoint(s) (either preselected or determined rigorously) ranked by
#'   significance and its/their associated significance value(s); Single value
#'   representing the MBIC threshold by which to base significance of a
#'   changepoint (i.e., penalty values LARGER than MBIC threshold indicate
#'   significance)
#'
#' @examples #Examples you wish to use, in demonstration of function
#' a = rpois(200,1)
#' b = rpois(200,2)
#' c = rpois(200,3)
#' v = c(a,b,c)
#'
#' change.poisson(v)
#' change.poisson(v, pre=400)
#'
#' @export
change.poisson = function(v, pre=NA) {
  if((sum(v<0)>0)){stop('Poisson test statistic requires positive data')}
  if(sum(as.integer(v)==v)!=length(v)){stop('Poisson test statistic requires integer data')}
  if(is.null(dim(v))==TRUE){
    n=length(v)
  }
  else{
    n=ncol(v)
  }
  if(n<4){stop('Data must have at least 4 observations to fit a changepoint model.')}
  if(n<(2*2)){stop('Minimum segment length is too large to include a change in this data')}
  if (is.na(pre) == FALSE) {
    cpt.pen = penalty_decision(penalty="MBIC", pen.value=0, n=n, diffparam=1, asymcheck="meanvar.poisson", method="AMOC")
    if(is.null(dim(v))==TRUE){
      tmp=preselect.poisson.calc(coredata(v),pre)
      tmp[3]=tmp[3]+log(tmp[1])+log(n-tmp[1]+1)
      my_row = list(c(tmp[1], tmp[2] - tmp[3]))
      my_df = as.data.frame(do.call(rbind,my_row))
      colnames(my_df) = c("Changepoint Location", "Changepoint Penalty Value")
      return(list(changepoint_data_frame = my_df, MBIC_penalty_threshold = cpt.pen))
    }
    else{
      tmp=preselect.poisson.calc(v,pre)
      tmp[,3]=tmp[,3]+log(tmp[,1])+log(n-tmp[,1]+1)
      rep=nrow(v)
      my_row=list()
      for(i in 1:rep){
        my_row[[i]]== c(i, tmp[i,1], tmp[i,2] - tmp[i,3])
      }
      my_df = as.data.frame(do.call(rbind,my_row))
      colnames(my_df) = c("Changepoint #", "Changepoint Location", "Changepoint Penalty Value")
      return(list(changepoint_data_frame = my_df, MBIC_penalty_threshold = cpt.pen))
    }
  } else {
    cpt = cpt.meanvar(v, test.stat="Poisson", method="BinSeg", minseglen=2, Q= n / 2 + 1)
    my_df = list()
    j = 1
    for (i in 1:(n/2 + 1)) {
      if (cpt@pen.value.full[i] >= cpt@pen.value) {
        my_row = c(i, cpt@cpts.full[i,i], cpt@pen.value.full[i])
        my_df[[j]] = my_row
        j = j + 1
      }
    }
    if (length(my_df) != 0) {
      my_df = as.data.frame(do.call(rbind,my_df))
      colnames(my_df) = c("Changepoint #", "Changepoint Location", "Changepoint Penalty Value")
      return(list(changepoint_data_frame = my_df, MBIC_penalty_threshold = cpt@pen.value))
    } else {
      return("No significant changepoints were found. If you wish to see the significance of a particular changepoint, please use the 'pre' parameter.")
    }
  }
}
