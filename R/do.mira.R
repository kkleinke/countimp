do.mira=function(imp,DV,fixedeff,randeff,grp,id,fam="poisson"){
  if (fam=="truncnbinom"){
    
    
    fit_zero<-vector(length=imp$m,"list")
    fit_count<-vector(length=imp$m,"list")
	  fe.copy<-fixedeff
      re.copy<-randeff


    
    for (j in 1:imp$m){
      data<-complete(imp,j)
      data[,id]<-as.factor(data[,id])
      data[,grp]<-as.factor(data[,grp])
      data$nz <- as.numeric(data[,DV]==0)
      
      if (length(fe.copy)>1){fi<-fe.copy[1]}else{fi<-fe.copy}
      if (length(re.copy)>1){ra<-re.copy[1]}else{ra<-re.copy}
      fit_zero[[j]] <- try(glmmPQL(fixed=as.formula(paste("nz~",fi,sep="")),random=as.formula(paste("~",ra,"|",grp,sep="")),
                                   data=data,
                                   family="binomial"))

      if (length(fe.copy)>1){fi<-fe.copy[2]}else{fi<-fe.copy}
      if (length(re.copy)>1){ra<-re.copy[2]}else{ra<-re.copy}

      fit_count[[j]] <-try(glmmadmb(formula=as.formula(paste(DV,"~",fi,sep="")),random=as.formula(paste("~",ra,"|",grp,sep="")),
                                    data=subset(data,data[,DV]>0),
                                    family="truncnbinom")) 
    }
est_fe.zero <- vector(length=imp$m,"list")
se_fe.zero <- vector(length=imp$m,"list")
est_fe.count <- vector(length=imp$m,"list")
se_fe.count <- vector(length=imp$m,"list")
est_re.zero <- vector(length=imp$m,"list")
est_re.count <- vector(length=imp$m,"list")
Rcor=vector(imp$m,mode="list")

for (j in 1:imp$m){  
  if (all(class(fit_zero)!="try-error")){
    est_fe.zero[[j]]<-fixef(fit_zero[[j]])
    se_fe.zero[[j]]<-summary(fit_zero[[j]])$tTable[,2]
    est_fe.count[[j]]<-fixef(fit_count[[j]])
    se_fe.count[[j]]<-summary(fit_count[[j]])$stdbeta
    est_fe.count[[j]]<-c(est_fe.count[[j]],alpha=fit_count[[j]]$alpha)  	
    se_fe.count[[j]]<-c(se_fe.count[[j]],alpha=fit_count[[j]]$sd_alpha)		
    est_re.count[[j]] <- sqrt(diag(VarCorr(fit_count[[j]])[[grp]])) #sd
    nam<-rownames(VarCorr(fit_zero[[j]]))
    est_re.zero[[j]] <- as.numeric(VarCorr(fit_zero[[j]])[,2])
    attr(est_re.zero[[j]],"names")<-nam
  }}
  miinf.zero<-data.frame(miinference(est_fe.zero, se_fe.zero))
  row.names(miinf.zero)<-paste(row.names(miinf.zero),".zero",sep="")
  miinf.count<-data.frame(miinference(est_fe.count, se_fe.count))
  row.names(miinf.count)<-paste(row.names(miinf.count),".count",sep="")
  MI<-rbind(miinf.zero,miinf.count)

  nam=names(est_re.count[[1]])
  R.count <- colMeans(matrix(unlist(est_re.count),nrow=imp$m,byrow=TRUE))
  attr(R.count,"names")<-paste(nam,".count",sep="")
  nam=names(est_re.zero[[1]])
  R.zero <- colMeans(matrix(unlist(est_re.zero),nrow=imp$m,byrow=TRUE))
  attr(R.zero,"names")<-paste(nam,".zero",sep="")
      MIRSD<-c(R.zero,R.count)
      
for (j in 1:imp$m){
	
	fit=fit_zero[[j]]
    fitsum=summary(fit)
    VC=VarCorr(fit)
	Rcor[[j]]=VC[,c(-1,-2)]
	    if(length(Rcor[[j]])>3){
      Rcor[[j]]=Rcor[[j]][-nrow(Rcor[[j]]),]
      coln=rownames(Rcor[[j]])
      
      suppressWarnings(mode(Rcor[[j]])<-"double")
      Rcor[[j]]=cbind(Rcor[[j]],NA)
      colnames(Rcor[[j]])<-coln
      Rcor[[j]][upper.tri(Rcor[[j]])]=t(Rcor[[j]][lower.tri(Rcor[[j]])])
      diag(Rcor[[j]])=1
    }else{
      if(length(Rcor[[j]])==3){
        Rcor[[j]]=Rcor[[j]][-length(Rcor[[j]])]
        coln=names(Rcor[[j]])
        
        suppressWarnings(mode(Rcor[[j]])<-"double")
        Rcor[[j]]=cbind(Rcor[[j]],NA)
        colnames(Rcor[[j]])<-coln
        Rcor[[j]][upper.tri(Rcor[[j]])]=t(Rcor[[j]][lower.tri(Rcor[[j]])])
        diag(Rcor[[j]])=1}else{Rcor[[j]]<-0}
    }
    if(any(Rcor[[j]]!=0)){
    colnames(Rcor[[j]])<-paste(colnames(Rcor[[j]]),".zero",sep="")      
rownames(Rcor[[j]])<-paste(rownames(Rcor[[j]]),".zero",sep="")      
 }
}
   
    
  MIRcor=Rcor[[1]]
  for (k in 2:imp$m)
  {
    MIRcor=MIRcor+Rcor[[k]]
  }
  
  MIRcor=MIRcor/imp$m

          output=list("FIT"=list("fit_zero"=fit_zero,"fit_count"=fit_count), "COEF.RA"=list("est_fe.zero"=est_fe.zero, "est_fe.count"=est_fe.count),
                "SE.RA"=list("se_fe.zero"=se_fe.zero, "se_fe.count"=se_fe.count),"MIINFERENCE"=MI,
                "R.SD"=list("est_re.zero"=est_re.zero, "est_re.count"=est_re.count),
                "R.COR"=Rcor, 
                "PR.SD"=MIRSD
                ,"PR.COR"=MIRcor
                )
    class(output) <- "do.mira"
    return(output)
  }else{
  
  
  # Declare temporary variables
  #############################
  
  m<-imp$m
  
  COEF=vector(imp$m,mode="list")	# fixed effects coefficients
  SE=vector(imp$m,mode="list")	# fixed effects SEs
  RSD=vector(imp$m,mode="list")	# random effets SDs
  Rcor=vector(imp$m,mode="list")	# random effects correlations
  FIT=vector(imp$m,mode="list")  # fitted glmmPQL model
  
  for (j in 1:imp$m)
  {
    cat("Imputation: ",j,"\n")
    tmp=complete(imp,j) # extract imputation # j from imp
    
    if(fam=="poisson"){
    fit <- glmmPQL(fixed=as.formula(paste(DV,"~",fixedeff,sep="")), 	data=tmp, 
                   random=as.formula(paste("~",randeff,"|",grp,sep="")),
                   family=fam,
                   control=list(opt="optim"), na.action=na.omit)
    FIT[[j]]<-fit}
    if(fam=="quasipoisson"){
      fit <- glmmPQL(fixed=as.formula(paste(DV,"~",fixedeff,sep="")),   data=tmp, 
                     random=as.formula(paste("~",randeff,"|",grp,sep="")),
                     family=fam,
                     control=list(opt="optim"), na.action=na.omit)
      FIT[[j]]<-fit}
    if(fam=="nbinom"){
      tmp[,id]<-as.factor(tmp[,id])
      tmp[,grp]<-as.factor(tmp[,grp])    
      fit0 <- glmmadmb(formula=as.formula(paste(DV,"~",fixedeff,sep="")),random=as.formula(paste("~",randeff,"|",grp,sep="")),
                       data=tmp,
                       family=fam)
      fit <- glmmPQL(fixed=as.formula(paste(DV,"~",fixedeff,sep="")),random=as.formula(paste("~",randeff,"|",grp,sep="")),
                     data=tmp,
                     family=negative.binomial(fit0$alpha),
                     control=list(opt="optim"))
      FIT[[j]]<-list(fit0,fit)}
      
    fitsum=summary(fit)
    VC=VarCorr(fit)
    
    COEF[[j]]=fitsum$tTable[,1]	# extract fixed effetcs: VALUE
    SE[[j]]=fitsum$tTable[,2]	# extract fixed effects: SE
    RSD[[j]]=VC[,2]
    mode(RSD[[j]])="double"
    Rcor[[j]]=VC[,c(-1,-2)]
    
    if(length(Rcor[[j]])>3){
      Rcor[[j]]=Rcor[[j]][-nrow(Rcor[[j]]),]
      coln=rownames(Rcor[[j]])
      
      suppressWarnings(mode(Rcor[[j]])<-"double")
      Rcor[[j]]=cbind(Rcor[[j]],NA)
      colnames(Rcor[[j]])<-coln
      Rcor[[j]][upper.tri(Rcor[[j]])]=t(Rcor[[j]][lower.tri(Rcor[[j]])])
      diag(Rcor[[j]])=1
    }else{
      if(length(Rcor[[j]])==3){
        Rcor[[j]]=Rcor[[j]][-length(Rcor[[j]])]
        coln=names(Rcor[[j]])
        
        suppressWarnings(mode(Rcor[[j]])<-"double")
        Rcor[[j]]=cbind(Rcor[[j]],NA)
        colnames(Rcor[[j]])<-coln
        Rcor[[j]][upper.tri(Rcor[[j]])]=t(Rcor[[j]][lower.tri(Rcor[[j]])])
        diag(Rcor[[j]])=1}else{Rcor[[j]]<-0}
    }
    
  } # end for (j in 1:imp$m)
   
  MI=miinference(COEF,SE)
  
  
  MIRSD=colMeans(t(as.data.frame(RSD)))
  
  MIRcor=Rcor[[1]]
  for (k in 2:imp$m)
  {
    MIRcor=MIRcor+Rcor[[k]]
  }
  
  MIRcor=MIRcor/imp$m
  
  if (fam=="nbinom"){
  theta <- vector(imp$m, mode="list")
  se.theta <- vector(imp$m, mode="list")
  for (i in 1:imp$m){theta[[i]]<-summary(FIT[[i]][[1]])$alpha}
  for (i in 1:imp$m){se.theta[[i]]<-summary(FIT[[i]][[1]])$sd_alpha}
  MI<-rbind(data.frame(MI),alpha=miinference(theta,se.theta))
  }
  output=list("FIT"=FIT, "COEF.RA"=COEF,"SE.RA"=SE,"MIINFERENCE"=data.frame(MI),"R.SD"=RSD,"R.COR"=Rcor, "PR.SD"=MIRSD,"PR.COR"=MIRcor)
  
  class(output) <- "do.mira"
  return(output)
  }#end else  
} # end function do.mira
