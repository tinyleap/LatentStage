#' Consideration Set Model
#'
#' Constructs consideration set via binary logit model, then selects from
#' consideration set via multinomial logit model.
#'
#' @param YObservedInput A vector of observed choices acting as response
#'   variable
#' @param XConsidInput A vector or matrix of choices potential to be included in
#'   consideration set
#' @param XSelectInput A vector or matrix of choices potential to be selected
#'   from consideration set (i.e., must be among \code{XConsidInput})
#' @param NumOpts Number of options for consideration set
#' @param ASCConsid Binary flag for whether users want ASC (i.e., intercept) in
#'   Consideration (i.e., choice). Default is 0.
#' @param ASCSelect Binary flag for whether users want ASC (i.e., intercept) in
#'   Selection (i.e., choice). Default is 0.
#' @param IncludedInput A vector or matrix of which choices are available;
#'   default is that every choice is available all the time.
#' @param meth Method of optimal function, default is method "L-BFGS-B". Other
#'   methods available are "BFGS", "CG", "SANN", "Brent".
#' @param Tolerance Controls the convergence of the "L-BFGS-B" method.
#'   Convergence occurs when the reduction in the objective is within this
#'   factor of the machine tolerance. Default is a tolerance of about 1e-5.
#'
#' @return A list object containing likelihood for consideration set model,
#'   likelihood for logit model, and the coefficient for all coviriates in
#'   Consideration and Selection. If we set ASCConsid/ASCSelect to be 1, the
#'   output will include ASC in Consideration/Selection.
#'
#' @examples
#' \donttest{
#' data(CFS5-1)
#' attach(CFS5-1)
#' RunConsiderationSet(V1, cbind(V2,V3,V4), cbind(V2,V3,V4), 5, 1, 1,Tolerance=1e-6)
#'
#' data(CFS6-1)
#' attach(CFS6-1)
#' RunConsiderationSet(V1, cbind(V2,V3,V4), cbind(V2,V3,V4), 6, 1, 1)
#'
#' data(CFS5-2)
#' attach(CFS5-2)
#' RunConsiderationSet(V1, cbind(V2,V3,V4), cbind(V2,V3,V4), 5, 1, 1, IncludedInput=V5, meth="BFGS")
#'
#' data(CFS6-2)
#' attach(CFS6-2)
#' RunConsiderationSet(V1, cbind(V2,V3,V4), cbind(V2,V3,V4), 6, 1, 1, IncludedInput=V5, meth="BFGS")
#' }
#'
#' @export RunConsiderationSet

RunConsiderationSet<-function(YObservedInput, XConsidInput, XSelectInput, NumOpts, ASCConsid=0, ASCSelect=0, IncludedInput=0*YObservedInput+1, meth="L-BFGS-B", Tolerance=1e-5){

  #if(nargs()<8)meth="L-BFGS-B"
  #if(nargs()<7)IncludedInput=0*YObservedInput+1
  #if(nargs()<6)ASCSelect=0
  #if(nargs()<5)ASCConsid=0
  Tolerance=Tolerance*1e15
  NumOccasions=length(YObservedInput)/NumOpts;

  # Make sure it works when there are only 1 variable input
  XConsidInput=as.matrix(XConsidInput)
  XSelectInput=as.matrix(XSelectInput)
  IncludedInput=as.matrix(IncludedInput)

  YObservedInputINIT=YObservedInput;
  XConsidInputINIT=XConsidInput;
  XSelectInputINIT=XSelectInput;
  IncludedInputINIT=as.matrix(IncludedInput);

  Y=t(matrix(YObservedInput, NumOpts, NumOccasions));
  YNew=NULL; XCNew=NULL; XSNew=NULL; IncNew=NULL; KeepInd=NULL;
  for(i in 1:NumOpts){
    A=Y; A[,i]=0;
    YNew=rbind(YNew, Y-A);
    XCNew=rbind(XCNew, XConsidInput);
    XSNew=rbind(XSNew, XSelectInput);
    IncNew=rbind(IncNew, IncludedInput);
    KeepInd=rbind(KeepInd, as.matrix(Y[,i]>0));
  }

  KeepMat=KeepInd%*%matrix(1,1,NumOpts);
  KeepIndFinal=t(KeepMat); KeepIndFinal=as.logical(as.numeric(KeepIndFinal));

  YNew=t(YNew); YNew=as.numeric(YNew);
  YObservedInput=YNew[KeepIndFinal];
  XConsidInput=XCNew[KeepIndFinal,];
  XSelectInput=XSNew[KeepIndFinal,];
  IncludedInput=IncNew[KeepIndFinal];

  XConsidInput=as.matrix(XConsidInput)
  XSelectInput=as.matrix(XSelectInput)

  NumOccasions=length(YObservedInput)/NumOpts;
  YObserved=t(matrix(YObservedInput, NumOpts, NumOccasions));
  TotalChoices=sum(YObservedInput)

  Included=t(matrix(IncludedInput, NumOpts, NumOccasions));


  # Test to be sure all chosen items are included in the set
  NotIncludedTest=sum(YObservedInput*(1-IncludedInput));
  # print(NotIncludedTest)
  if(NotIncludedTest>0){
    cat('Chosen Items Must Be Among Included Set');cat('\n')
    cat(NotIncludedTest);cat('\n')
    #return(NULL)
  }

  if(ASCConsid==1) XConsid=cbind(kronecker(matrix(1,NumOccasions,1), diag(NumOpts)),XConsidInput) else XConsid=XConsidInput

  if(ASCSelect==1) {
    XSelectASCTemp=kronecker(matrix(1,NumOccasions,1), diag(NumOpts));
    XSelect=cbind(as.matrix(XSelectASCTemp[,1:(NumOpts-1)]),XSelectInput);
  }else XSelect=XSelectInput;
  LengthXConsid = dim(XConsid)[1];
  WidthXConsid = dim(XConsid)[2];
  LengthXSelect=dim(XSelect)[1];
  WidthXSelect = dim(XSelect)[2];

  # This is for the ordinary Logit model
  res1 = MaxLogitLL(matrix(0,WidthXConsid+WidthXSelect,1), YObserved, XConsid, XSelect, meth, Included, Tolerance);
  BetaVectLogit=res1$BetaVectLogit
  fvalLogit=res1$fvalLogit
  exitflagLogit=res1$exitflagLogit
  gradLogit=res1$gradLogit
  hessianLogit=res1$hessianLogit

  # This is for the consideration set model
  B=as.matrix(c(0,1));
  A=B; for(i in 1:(NumOpts-1)) A=cbind(kronecker(as.matrix(c(0,1)), 1+as.matrix(0*A[,1])),kronecker(as.matrix(c(1,1)), A))
  AllSets=A;
  remove(A)
  remove(B)

  res2 = MaxConsidLL(matrix(0,WidthXConsid+WidthXSelect,1), YObserved, XConsid, XSelect, meth, AllSets, Included, Tolerance);
  BetaVect=res2$BetaVect
  fval=res2$fval
  exitflag=res2$exitflag
  grad=res2$grad
  hessian=res2$hessian

  LL=-fval;

  cat('\n')
  cat('TotalChoices = ')
  cat(TotalChoices)
  cat('\n')
  cat('\n')
  cat(paste('LogLikelihood, Consideration + Selection = ',as.character(round(LL,4)),'\n'))
  cat(paste('LogLikelihood, Selection Only            = ',as.character(round(-fvalLogit,4)),'\n'))
  cat(paste('Chi-Square Likelihood Ratio Test Value   = ',as.character(round(2*(fvalLogit-fval),4)) ,', with ' ,as.character(WidthXConsid),' df','\n'))
  cat(paste('p-value on excluding Consideration       = ',sprintf("%.4e",round(1-pchisq(2*(fvalLogit-fval),WidthXConsid),20)),'\n'))
  cat(paste('% Improvement due to CS per Choice       = ',as.character(round(100*(exp((fvalLogit-fval)/TotalChoices)-1),2)),' %','\n'))
  cat('\n')

  StdErrs=sqrt(diag(solve(hessian)));
  StdErrs=as.matrix(StdErrs,length(StdErrs),1)
  ZValues = abs(BetaVect/StdErrs);
  PValues = 2*pnorm(-abs(ZValues));
  OutputMatrix=cbind(BetaVect,StdErrs,ZValues,PValues);

  if(ASCConsid==1){cat('ASCs, Consideration\n')
    mat = round(as.matrix(OutputMatrix[1:NumOpts,]),4)
    if(NumOpts==1){
      p(t(mat))
    } else p(mat)

    cat(' \n')
    cat('Coeffs, Consideration\n')
    mat = round(as.matrix(OutputMatrix[(NumOpts+1):(WidthXConsid),]),4)
    if((WidthXConsid-NumOpts)==1){
      p(t(mat))
    } else p(mat)
    cat('\n')
  }

  if (ASCConsid==0){
    cat('Coeffs, Consideration\n')
    mat = round(as.matrix(OutputMatrix[1:WidthXConsid,]),4)
    if(WidthXConsid==1){
      p(t(mat))
    } else p(mat)
    cat('\n')
  }

  if (ASCSelect==1){
    cat('ASCs, Selection\n')
    mat = round(as.matrix(OutputMatrix[(WidthXConsid+1):(WidthXConsid+NumOpts-1),]),4)
    if(NumOpts==2){
      p(t(mat))
    } else p(mat)

    cat(' \n')
    cat('Coeffs, Selection\n')
    mat = round(as.matrix(OutputMatrix[(WidthXConsid+NumOpts):(WidthXConsid+WidthXSelect),]),4)
    if(WidthXSelect==NumOpts){
      p(t(mat))
    } else p(mat)
    cat(' \n')
  }

  if (ASCSelect==0){
    cat('Coeffs, Selection\n')
    mat = round(as.matrix(OutputMatrix[(WidthXConsid+1):(WidthXConsid+WidthXSelect),]),4)
    if(WidthXSelect==1){
      p(t(mat))
    } else p(mat)
    cat(' \n')
  }
  cat(' \n')
  # return(list=list(LL=LL, BetaVect=BetaVect, StdErrs=StdErrs, ZValues=ZValues, PValues=PValues, OutputMatrix=OutputMatrix, Hessian=hessian))
}

MaxLogitLL<-function(BetaBoth, YObserved, XConsid, XSelect, meth, Included, Tolerance){
  NumOccasions=dim(YObserved)[1];
  NumOpts=dim(YObserved)[2]
  if(nargs()<5)Included=matrix(1,NumOccasions,NumOpts)
  NotIncludedTest=sum(YObserved*(1-Included));
  #if(NotIncludedTest>0) return(NULL);


  XConsidLength=dim(XConsid)[1];
  XConsidWidth=dim(XConsid)[2]
  XSelectLength=dim(XSelect)[1];
  XSelectWidth=dim(XSelect)[2]

  LogitLL1<-function(x){
    BetaSelect=x[(XConsidWidth+1):(XConsidWidth+XSelectWidth)];
    BetaSelect=matrix(BetaSelect,XSelectWidth,1)
    USelect=XSelect%*%BetaSelect;
    PSelectRaw=exp(USelect);
    PSelectNewRaw=t(matrix(PSelectRaw,NumOpts,NumOccasions));
    IncludedNew=t(matrix(t(Included),NumOpts,NumOccasions));

    ChosenNewInit=t(matrix(t(YObserved),NumOpts,NumOccasions));
    Multiplicities=rowSums(ChosenNewInit);
    ChosenNew=(ChosenNewInit>0);

    YNew=-sum(log(colSums(t(PSelectNewRaw)*t(IncludedNew)*t(ChosenNew)))*t(Multiplicities))+sum(log(colSums(t(PSelectNewRaw)*t(IncludedNew)))*t(Multiplicities));
    #    Intermediate=cbind(colSums(IncludedNew), colSums(YObservedNew), colSums(PSelectNewRaw));
    #    print(as.matrix(Intermediate))
    Y=YNew;
    return(Y)
  }
  res2=optim(BetaBoth, LogitLL1, gr = NULL, method = meth,control = list(maxit=100000,factr=Tolerance), hessian = TRUE)
  Z=res2$par
  fval=res2$value
  grad=res2$counts[2]
  hessian=res2$hessian
  exitflag=res2$convergence
  return(list=list(BetaVectLogit=Z,fvalLogit=fval,exitflagLogit=exitflag,gradLogit=grad,hessianLogit=hessian))
}


MaxConsidLL<-function(BetaBoth, YObserved, XConsid, XSelect, meth, AllSets, Included, Tolerance){
  NumOccasions=dim(YObserved)[1];
  NumOpts=dim(YObserved)[2]
  if(nargs()<5)Included=matrix(1,NumOccasions,NumOpts)
  NotIncludedTest=sum(YObserved*(1-Included));
  #if(NotIncludedTest>0) return(NULL);

  XConsidLength=dim(XConsid)[1];
  XConsidWidth=dim(XConsid)[2]
  XSelectLength=dim(XSelect)[1];
  XSelectWidth=dim(XSelect)[2]

  ConsidLL1<-function(x){
    BetaConsid=x[1:XConsidWidth];
    BetaSelect=x[(XConsidWidth+1):(XConsidWidth+XSelectWidth)];

    UConsid=XConsid%*%matrix(BetaConsid,XConsidWidth,1);
    # PConsidNew is a matrix of the Consideration probabilities for each option (columns) for each choice occasion (rows)
    PConsid=1/(1+exp(-UConsid));
    PConsidNew=t(matrix(PConsid,NumOpts,NumOccasions));

    USelect=XSelect%*%matrix(BetaSelect,XSelectWidth,1);
    USelectNew=t(matrix(t(USelect),NumOpts,NumOccasions));
    USelectNewRescaled=USelectNew-matrix(apply(t(USelectNew),2,max),ncol=1)%*%matrix(1,1,NumOpts); # Subtracts largest utility in each row, helps stability
    expUSelectNewRescaled=exp(USelectNewRescaled);

    ChosenNewInit=t(matrix(t(YObserved),NumOpts,NumOccasions));
    IncludedNew=t(matrix(t(Included),NumOpts,NumOccasions));
    Multiplicities=t(colSums(t(ChosenNewInit))); # This stores HOW MANY they chose
    ChosenNew=(ChosenNewInit>0); # This stores WHAT they chose, not how many

    A=PConsidNew; # Consideration Set probabilities for each option
    B=expUSelectNewRescaled; # exp(Beta*X) for selection
    C=ChosenNew; # Items chosen
    D=IncludedNew; # Items in set
    Sets=t(AllSets);

    ANew=A*D + .00000000001*(1-D);
    # CSNum=exp(log(A)%*%Sets+log(1-A)%*%(1-Sets)); # Probability of each CS for each occasion
    CSNum=exp(log(ANew)%*%Sets+log(1-ANew)%*%(1-Sets)); # Probability of each CS for each occasion
    # print(as.matrix(rowSums(CSNum))); # checks that the rows sum to 1
    SelectDen=B%*%Sets;
    SelectNum=(B*C)%*%Sets;
    SelectAll=SelectNum/SelectDen;
    Final1=sum(log(colSums(t(CSNum[,-1])*t(SelectAll[,-1])))*t(Multiplicities));
    Final2=-sum(log(1-CSNum[,1])*Multiplicities);

    LLTotal=Final1+Final2;

    Y=-LLTotal;
    return(Y)
  }


  # res2=optim(BetaBoth, ConsidLL1, gr = NULL, method = "BFGS",control = list(maxit=10000, reltol=.000001), hessian = TRUE)
  # res2=optim(BetaBoth, ConsidLL1, gr = NULL, method = "CG",control = list(maxit=10000, reltol=.000001), hessian = TRUE)

  res2=optim(BetaBoth, ConsidLL1, gr = NULL, method = meth, control = list(maxit=100000,factr=Tolerance), hessian = TRUE)
  Z=res2$par
  fval=res2$value
  grad=res2$counts[2]
  hessian=res2$hessian
  exitflag=res2$convergence
  return(list=list(BetaVect=Z,fval=fval,exitflag=exitflag,grad=grad,hessian=hessian))
}
p <- function(x){
  prmatrix(x, rowlab=rep("",nrow(x)), collab=c("Coeff","StdErr","ZValue","PValue"))
}
