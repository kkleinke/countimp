mice.impute.2l.nb2.noint.boot <-
function(y,ry,x,type,intercept=FALSE)
{
  # allowed predictorMatrix entries:
  # class variable = -2 (only one allowed!)
  # random and fixed effect = 2
  # fixed effect only = 1
  # variable NOT in imputation model = 0
  
  require(glmmADMB)
  if (!require(pscl)){stop("glmmADMB package must be 
  installed!")}
  
  Y=y[ry]
  X=x[ry,]
  
  nam=paste("V",1:ncol(X),sep="")
  colnames(X)=nam
  if(sum(type==-2)>1){stop("only one class allowed!")}
  
  grp=which(type==-2)
  X[,grp]<-as.factor(X[,grp])
  
  fixedeff=paste("+", paste(nam[-grp],collapse="+"),sep="")
  if (any(type==2)){
    
    ran=which(type==2)
    
    randeff=paste("+",paste(nam[ran],collapse="+"),sep="")
    if(!intercept)
    {
      randeff=paste("~0",randeff,"|",paste(nam[grp]),sep="")
    } else {
      randeff=paste("~1",randeff,"|",paste(nam[grp]),sep="")
    }
  }else{
    if(!intercept)
    {
      randeff <- paste("~0","|",paste(nam[grp]),sep="")
    } else {
      randeff <- paste("~1","|",paste(nam[grp]),sep="")
    }  
  }
  randeff=as.formula(randeff)
  
  fixedeff=as.formula(paste("Y","~",fixedeff,sep=""))
  dat=data.frame(Y,X)
  datBS <- dat[sample(1:length(Y),length(Y),replace=TRUE),]
  fit<-  
    glmmadmb(
      formula=fixedeff,
      random=randeff,
      data=datBS,
      family="nbinom")
  
  library(MASS)
  fit2<- glmmPQL(fixed=fixedeff,random=randeff, data=datBS,
                 family=negative.binomial(fit$alpha))
  
  
  newdatamis=data.frame(X=x[!ry,])
  colnames(newdatamis)=nam	
  
  p<- predict(fit2,newdata=newdatamis, type="response", na.action=na.pass)
  
  
  im=rnegbin(length(p),theta=fit$alpha,mu=p) 
  
  return(im)
}
