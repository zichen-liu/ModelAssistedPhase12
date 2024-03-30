rm(list=ls())
ndose<-6 # number of doses
targetT<-0.3 # highest acceptable toxicity rate
targetE<-0.2 # lowest acceptable efficacy rate
u1<-100
u2<-0
ntrial<-10000
obd.vec<-mtd.vec<-c()
pE.mat<-pT.mat<-u.mat<-c()
pE.obd<-pT.obd<-c()
pE.mtd<-pT.mtd<-c()
set.seed(1)
for(s in 1:ntrial){
  
  obd<-sample(ndose,1)
  mtd<-obd # initialization
  prob_plateau<-0.5
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
  
  pE.mat<-rbind(pE.mat,jj)
  pT.mat<-rbind(pE.mat,kk)
  u.mat<-rbind(u.mat,utility)
  obd.vec<-c(obd.vec,obd)
  mtd.vec<-c(mtd.vec,mtd)
  
  pE.obd<-c(pE.obd,jj[obd])
  pT.obd<-c(pT.obd,kk[obd])
  
  pE.mtd<-c(pE.mtd,jj[mtd])
  pT.mtd<-c(pT.mtd,kk[mtd])

}

table(obd.vec)
hist(pE.obd)
quantile(pE.obd)
hist(pT.obd)
quantile(pT.obd)

