mice.impute.2l.zihnb.boot <-
function(y,ry,x,type,intercept.c=TRUE,intercept.z=TRUE){
  # type is extracted from predictorMatrix
  # y: variable to be imputed
  # x: predictors
  # ry: response indicator
  # y,ry,x are passed on by mice  
  
  library(MASS) # needed for glmmPQL model
  Y=y[ry]
  X=x[ry,]
  X=data.frame(X)	
  nam=colnames(X)
  # only two-level model allowed / one class variable!
  if(sum(type==-2)>1){stop("only one class allowed!")}
  grp<-which(type==-2) # group variable
  fc<-which(type==3) #fixed only count model
  rc<-which(type==4) #random and fixed count model
  fz<-which(type==5) #fixed only zero model
  rz<-which(type==6) #random and fixed zero model
  f<-which(type==1) #fixed only both models
  r<-which(type==2) #random and fixed both model
  X[,grp]<-as.factor(X[,grp]) # group variable must be categorical!
  
  #zero model
  ###########
  
   # formula for random part of the model
  ran<-c(r,rz); ran<-unique(ran);ran<-sort(ran)
  randeff<-paste(nam[ran],collapse="+")
 if (any(type==2)|any(type==6)){
  	if (intercept.z==TRUE){
  	randeff<-paste("~1+",randeff,"|",paste(nam[grp]),sep="")
  	}else{
  	    randeff<-	
  	    paste("~0+",randeff,"|",paste(nam[grp]),sep="")	
  		}
  }else{
  	if (intercept.z==TRUE){
  	randeff <- paste("~1","|",paste(nam[grp]),sep="")
  	}else{
  		#randeff <- paste("~0","|",paste(nam[grp]),sep="")
  		stop("at least one random effect needed")
  		}
  	}
 randeff<-as.formula(randeff)
  
  #formula for fixed part of the model
  fix<-c(f,r,fz,rz);fix<-unique(fix);fix<-sort(fix)
  fixedeff<-paste(nam[fix],collapse="+")
  fixedeff<-as.formula(paste("nz","~",fixedeff,sep=""))	
  nz<-as.numeric(Y==0)
  dat<-data.frame(Y,X,nz)	
  datBS <- dat[sample(1:length(Y),length(Y),replace=TRUE),]
  # zero model is a 2L binomial glmm
  fit <-
    glmmPQL(
      fixed=fixedeff,
      random=randeff,
      data=datBS,
      family="binomial")
  
  # Bayesian regression  
  
  newdatamis=data.frame(X=x[!ry,])
  colnames(newdatamis)=nam	
  
  #compute predicted probabilities for having 
  # zero vs. non-zero count
  p <- predict(fit,
               newdata=newdatamis,type="response",na.action=na.pass)
  # imputed values are drawn from binomial distribution
  im<-rbinom(n=length(p),size=1,prob=1-p)
  
  nonzero<-which(im==1) # which cases have non-zero values
  
  
  # count model
  #############
  
  library(glmmADMB) # needed to fit model
  
  # random effects formula
  ran<-c(r,rc);ran<-unique(ran);ran<-sort(ran)
  randeff<-paste(nam[ran],collapse="+")
 if (any(type==2)|any(type==4)){
  	if (intercept.c==TRUE){
  	randeff<-paste("~1+",randeff,"|",paste(nam[grp]),sep="")
  	}else{
  	    randeff<-	
  	    paste("~0+",randeff,"|",paste(nam[grp]),sep="")	
  		}
  }else{
  	if (intercept.c==TRUE){
  	randeff <- paste("~1","|",paste(nam[grp]),sep="")
  	}else{
  		#randeff <- paste("~0","|",paste(nam[grp]),sep="")
  		stop("at least one random effect needed")
  		}
  	} 
  randeff<-as.formula(randeff)
  
  # fixed effects formula	
  fix<-c(f,r,fc,rc);fix<-unique(fix);fix<-sort(fix)
  fixedeff<-paste(nam[fix],collapse="+")
  fixedeff<-as.formula(paste("Y","~",fixedeff,sep=""))	
  
  #count model is a zero-truncated 2L NB model
  fit<-	
    glmmadmb(
      formula=fixedeff,
      random=randeff,
      data=subset(datBS,datBS[,"Y"]>0),
      family="truncnbinom")
  
  # Bayesian regression    
  newdatamis=data.frame(X=x[!ry,])
  colnames(newdatamis)=nam	
  
  # predict function of glmmADMB package 
  # can only handle FIXED effects
  # and must thus NOT be used!
  
  fit.mat <- matrix(0, nrow(newdatamis), 2, 
                    dimnames=list(1:nrow(newdatamis),c("fixed","random")))
  
  nf=names(fixef(fit))
  nf=nf[-which(nf=="(Intercept)")]
  
  Xnew=data.matrix(cbind(1,newdatamis[,nf]))
  
  
  
  
  #predictions based on fixed effects only
  fit.mat[,"fixed"] <- Xnew%*%fixef(fit) 
  
  # restructure data set
  nr=colnames(ranef(fit)[[1]])
  nr=nr[-which(nr=="(Intercept)")]
  Xnew2=data.matrix(cbind(1,newdatamis[,nr]))
  RE=ranef(fit)[[1]]
  RE=data.matrix(RE[newdatamis[,grp],])
  
  # predictions based on fixed + random effects!
    
  fit.mat[,"random"] <- 
    fit.mat[, "fixed"] + data.matrix(Xnew2*RE)%*% rep(1,ncol(RE))
    
  library(aster)
  # imputed values drawn from zero-truncated NB distribution
  Yfit_p2=rktnb(nrow(fit.mat), 
                size=fit$alpha, k=0, mu=exp(fit.mat[,"random"]), xpred = 1) 
  im[nonzero]<-Yfit_p2[nonzero]
  return(im)
}
