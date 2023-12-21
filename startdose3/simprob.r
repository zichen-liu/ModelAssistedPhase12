
simprob<-function(ndose,targetE,targetT,u1,u2 ,randomtype){
  
  obd<-sample(ndose,1)
  mtd<-obd # initialization
  if(randomtype == 1){
  prob_plateau<-0.5
  }else{
  prob_plateau<-0
  }
  obd.temp<-0
  jj<-kk<-rep(0,ndose)
  uu<-runif(1)
  
  if(uu<prob_plateau & obd<ndose){
    # generate plateau cases 
    while(obd.temp!=obd | kk[obd]>(targetT+0.1) | jj[obd]<targetE){
      mtd<-(obd-1)+sample(ndose-obd+1,1)
      D=1:ndose
      if(mtd<ndose){temp<-sort(runif(length(D)-mtd,targetT,1));
      bornesup<-max(targetT+0.05,temp[length(D)-mtd])}
      if(mtd==ndose){bornesup<-targetT +(1-targetT)*rbeta(1,0.5,1)}
      mtd.temp<-0
      
      while(mtd.temp!=mtd | kk[mtd]< (targetT-0.05)){
        kk<-sort(runif(ndose,0,bornesup))
        mtd.temp<-which.min(abs(kk-targetT))
      }
      
      bornesup<-runif(1,0.3,0.8)
      jj[1:obd]<-sort(runif(obd,0,bornesup))
      jj[obd:ndose]<-jj[obd]
      
      utility<-u1*jj+(1-kk)*u2
      
      obd.temp<-which.max(utility[1:mtd])
      
    }
  } 
  
  
  if(uu>=prob_plateau & obd<ndose){
    # generate non-plateau cases 
    while(obd.temp!=obd | kk[obd]>(targetT+0.1) | jj[obd]<targetE){
      mtd<-(obd-1)+sample(ndose-obd+1,1)
      D=1:ndose
      if(mtd<ndose){temp<-sort(runif(length(D)-mtd,targetT,1));
      bornesup<-max(targetT+0.05,temp[length(D)-mtd])}
      if(mtd==ndose){bornesup<-targetT +(1-targetT)*rbeta(1,0.5,1)}
      mtd.temp<-0
      
      while(mtd.temp!=mtd | kk[mtd]< (targetT-0.05)){
        kk<-sort(runif(ndose,0,bornesup))
        mtd.temp<-which.min(abs(kk-targetT))
      }
      
      med<-(obd-1)+sample(ndose-obd+1,1)
      bornesup<-runif(1,0.3,0.8)
      jj[med]<-max(runif(ndose,0,bornesup))
      if(med>1){
        jj[1:(med-1)]<-sort(runif(med-1,0,jj[med]))
      }
      if(med<ndose){
        jj[(med+1):ndose]<- sort(runif(ndose-med,0,jj[med]),,decreasing=TRUE)
      }
      utility<-u1*jj+(1-kk)*u2
      obd.temp<-which.max(utility[1:mtd])
      
    }
  }
  
  if(obd==ndose){
    while(obd.temp!=obd | kk[obd]>(targetT+0.1) | jj[obd]<targetE){
      mtd<-(obd-1)+sample(ndose-obd+1,1)
      D=1:ndose
      if(mtd<ndose){temp<-sort(runif(length(D)-mtd,targetT,1));
      bornesup<-max(targetT+0.05,temp[length(D)-mtd])}
      if(mtd==ndose){bornesup<-targetT +(1-targetT)*rbeta(1,0.5,1)}
      mtd.temp<-0
      
      while(mtd.temp!=mtd | kk[mtd]< (targetT-0.05)){
        kk<-sort(runif(ndose,0,bornesup))
        mtd.temp<-which.min(abs(kk-targetT))
      }
      bornesup<-runif(1,0.3,0.8)
      jj<-sort(runif(ndose,0,bornesup))
      utility<-u1*jj+(1-kk)*u2
      obd.temp<-which.max(utility[1:mtd])
    }
  }
  
  return(list(pE = jj,pT=kk))
  
}
