setwd('/tmp/temp')


ndose <<- 5 
OBD <<- 0;

for (iii in 1:ndose){
OBD <<- iii 


###############################################################
###############################################################
##########   simulate nsimul tirals   #########################
ji3plus3_sim <- function(pT=0.3,  pE=0.2, eps1=0.05, eps2=0.05, 
                         psafe=0.95, pfutility=0.95, 
                         cohortsize=3,ncohort=10,startdose=1,
                         nsimul=1000,utilitytype = 1, randomtype = 1){
  
  
  if (utilitytype==1){
    u1 = 60 
    u2 = 40 
  }
  if (utilitytype==2){
    u1 = 100 
    u2 = 0 
  }
  
  
  
  #########################################
  # Obtain optimal dosing decsion for ji3+3
  ji3plus3 <- function(x, y, n, pT, eps1, eps2, pE){
    if ((x/n < pT-eps1) && (y/n <= pE)){
      dec <- "E"
    } else if ((x/n < pT-eps1) && (y/n > pE)){
      dec <- "S"
    } else if ((x/n >= pT-eps1 & x/n <= pT+eps2) && (y/n <= pE)){
      dec <- "E"
    } else if ((x/n >= pT-eps1 & x/n <= pT+eps2) && (y/n > pE)){
      dec <- "S"
    } else if ((x/n > pT+eps2) && ((x-1)/n < pT-eps1)) {
      dec <- "S"
    } else {
      dec <- "D"
    }
    return(dec)
  }
  
  
  library(Iso)
  
  
  simprob<-function(ndose,targetE,targetT,u1,u2 ,randomtype){
    targetE=0.20
    targetT=0.30
    
    if(OBD==0){
      obd<-sample(ndose,1)
    }else{
      obd <- OBD
    }
    
    mtd<-obd # initialization
    prob_plateau<-0.5
    obd.temp<-0
    jj<-kk<-rep(0,ndose)
    uu<-runif(1)
    
    if(uu<prob_plateau & obd<ndose){
      # generate plateau cases 
      while(obd.temp!=obd | kk[obd]>(targetT+0.05) | jj[obd]<(targetE+0.05)){
        mtd<-(obd-1)+sample(ndose-obd+1,1)
        D=1:ndose
        if(mtd<ndose){temp<-sort(runif(length(D)-mtd,targetT,1));
        bornesup<-max(targetT+0.10,temp[length(D)-mtd])}
        if(mtd==ndose){bornesup<-targetT +(1-targetT)*rbeta(1,0.5,1)}
        mtd.temp<-0
        
        while(mtd.temp!=mtd){
          kk<-sort(runif(ndose,0,bornesup))
          mtd.temp<-which.min(abs(kk-targetT))
        }
        
        bornesup<-runif(1,targetE+0.10,1.0)
        #med=obd
        #bornesup<-max(runif(med,targetE+0.10,0.9))
        jj[1:obd]<-sort(runif(obd,0,bornesup))
        jj[obd:ndose]<-jj[obd]
        
        utility<-u1*jj+(1-kk)*u2
        
        obd.temp<-which.max(utility[1:mtd])
        
      }
    } 
    
    
    if(uu>=prob_plateau & obd<ndose){
      # generate non-plateau cases 
      while(obd.temp!=obd | kk[obd]>(targetT+0.05) | jj[obd]<(targetE+0.05)){
        mtd<-(obd-1)+sample(ndose-obd+1,1)
        D=1:ndose
        if(mtd<ndose){temp<-sort(runif(length(D)-mtd,targetT,1));
        bornesup<-max(targetT+0.10,temp[length(D)-mtd])}
        if(mtd==ndose){bornesup<-targetT +(1-targetT)*rbeta(1,0.5,1)}
        mtd.temp<-0
        
        while(mtd.temp!=mtd){
          kk<-sort(runif(ndose,0,bornesup))
          mtd.temp<-which.min(abs(kk-targetT))
        }
        
        med<-(obd-1)+sample(ndose-obd+1,1)
        bornesup<-runif(1,targetE+0.10,1.0)
        #bornesup<-max(runif(med,targetE+0.10,0.9))
        jj[med]<-max(runif(ndose,0,bornesup))
        if(med>1){
          jj[1:(med-1)]<-sort(runif(med-1,0,jj[med]))
        }
        if(med<ndose){
          jj[(med+1):ndose]<- sort(runif(ndose-med,0,jj[med]),decreasing=TRUE)
        }
        utility<-u1*jj+(1-kk)*u2
        obd.temp<-which.max(utility[1:mtd])
        
      }
    }
    
    if(obd==ndose){
      while(obd.temp!=obd | kk[obd]>(targetT+0.05) | jj[obd]<(targetE+0.05)){
        mtd<-(obd-1)+sample(ndose-obd+1,1)
        D=1:ndose
        if(mtd<ndose){temp<-sort(runif(length(D)-mtd,targetT,1));
        bornesup<-max(targetT+0.10,temp[length(D)-mtd])}
        if(mtd==ndose){bornesup<-targetT +(1-targetT)*rbeta(1,0.5,1)}
        mtd.temp<-0
        
        while(mtd.temp!=mtd){
          kk<-sort(runif(ndose,0,bornesup))
          mtd.temp<-which.min(abs(kk-targetT))
        }
        bornesup<-runif(1,targetE+0.10,1.0)
        #med=ndose
        #bornesup<-max(runif(med,targetE+0.10,0.9))
        jj<-sort(runif(ndose,0,bornesup))
        utility<-u1*jj+(1-kk)*u2
        obd.temp<-which.max(utility[1:mtd])
      }
    }
    
    return(list(pE = jj,pT=kk,obd=obd,mtd=mtd))
  }
  
  
  
  peestimate<-function(yE,n){
    ndose<-length(yE)
    lik<-rep(0,ndose)
    pe<-(yE+0.05)/(n+0.1)
    p.e<-matrix(NA,ndose,ndose)
    for (i in 1:ndose){
      if (i==1) {x<-seq(ndose,1,by=-1)} else {x<-c(1:(i-1),seq(ndose,i))}
      #x<-x
      p.e[i,]<-ufit(pe,lmode=i,x=x,w=n+0.5)[[2]]
      lik[i]<-prod(dbinom(yE,n,p.e[i,]))		
    }
    lik<-lik/sum(lik)
    pe<-t(p.e)%*%lik+0.01*seq(1,ndose)
    return(pe)}
  
  
  
  #########################################
  ## pava is the pool adjacent violator algorithm to perform isotonic 
  ## transformation to get the ordered posterior sample
  pava <- function (x, wt = rep(1, length(x))) 
  {
    n <- length(x)
    if (n <= 1) 
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol))) 
        break
      i <- min((1:(n - 1))[viol])
      # print(i)
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      # print(ilvl)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      # print(x)
      lvlsets[ilvl] <- lvl1
      # print(lvlsets)
    }
    x
  }
  
  #########################################
  ## betavar computes variances of beta distributions 
  betavar<-function(a,b){
    resp <- a*b/((a+b)^2*(a+b+1))
    return(resp)
  }
  
  
  
  
  qE <- pE
  
  parabeta <- c(1,1) #Beta distribution parameters
  
  selection <- numeric(nsimul)
  pat_treated <- matrix(NA, nrow = nsimul, ncol = ndose)
  
  
  csize= cohortsize
  samplesize = csize * ncohort 
  
  targetT = pT
  targetE = 0.20
  ntrial = nsimul
  npts = ncohort*cohortsize;
  
  
  bd.sel=0;
  bd.pts=0;
  od.sel=0;
  od.pts=0;
  ov.sel=0;
  ntox=0;
  neff=0;
  poorall=0;
  incoherent=0;
  overdose=0;
  
  
  dselect = rep(0, ntrial); # store the selected dose level
  
  set.seed(30);
  for (simul in 1:nsimul){ 
    
    
    probs = simprob(ndose,targetE,targetT,u1,u2 ,randomtype)
    jj = probs$pE
    kk = probs$pT
    pE.true<-jj
    pT.true<-kk
    u.true<-(u1*pE.true+(1-pT.true)*u2);
    #bd=which.max(u.true*(pT.true<(targetT+0.1))*(pE.true>(targetE-0.05)));
    bd = probs$obd
    mtd = probs$mtd
    toxtrue <- pT.true
    efftrue <- pE.true
    
    
    x <- rep(0,ndose) #toxicity number for each dose
    y <- rep(0,ndose) #efficacy number for each dose
    n <- rep(0,ndose) #patient number for each dose
    d <- startdose #current dose is the lowest dose at the beginning of the trial
    st <- 0 # stop sign
    currentsize <- 0
    dtox <- rep(0,ndose)
    dfut <- rep(0,ndose)
    
    earlystop = 0 
    
    while(st == 0){
      xx <-sum(runif(cohortsize)<toxtrue[d]) # generate tox outcome at current dose
      yy <-  sum(runif(cohortsize)<efftrue[d]) # generate eff outcome at current dose
      x[d] <- x[d] + xx # total tox outcome
      y[d] <- y[d] + yy # total eff outcome
      n[d] <- n[d] + csize # total sample size used
      currentsize <- currentsize + csize
      # safety rule	(monotonic tox)	
      if (pbeta(pT, parabeta[1] + x[d], parabeta[2] + n[d] - x[d]) < 1 - psafe) {
        dtox[d:ndose] <- 1
      }
      # futility rule (mark only current dose)
      if (pbeta(0.20, parabeta[1] + y[d], parabeta[2] + n[d] - y[d]) > 0.90) {
        dfut[d] <- 2
      }
      # obtain optimal dosing decisions
      dec <- ji3plus3(x=x[d], y=y[d], n=n[d], pT=pT, eps1=eps1, eps2=eps2, pE=pE)
      # if the sampsize is reached, stop
      if (currentsize >= samplesize){
        st <- 1
      }
      else {
        # if all doses are either too toxic or of no efficacy (couldn't find a dose for the next corhort)
        if (length(which(dfut+dtox==0))==0) {
          st <- 1 
          earlystop = 1
        }
        else {
          # if current dose is too toxic
          if (dtox[d]==1){
            # if there is a lower dose available, de-escalate to the closet dose below
            if (d!=1 && length(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0)!=0)) d <- max(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0)) 
            else {
              st <- 1 # if there is no valid dose below current dose
            }
          }
          # if current dose is of no efficacy
          else if (dfut[d]==2) { 
            if (dec=="E") { 
              # if there is a valid dose above current dose, escalate by minimum dose size
              if (d!=ndose && length(which(dtox[(d+1):ndose]+dfut[(d+1):ndose]==0))!=0) d <- min(which(dtox[(d+1):ndose]+dfut[(d+1):ndose]==0))+d 
              # if there is no dose above current dose, de-escalate to the next available dose
              else if (d!=1 && length(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0)!=0)) d <- max(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0)) 
              # if there is no dose above or below current dose, terminate the trial
              else {
                st <- 1
                earlystop = 1
              }
            }
            else if (dec=="D") {
              # if there is a valid dose below current dose, de-escalate by minimum dose size
              if (d!=1 && length(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0)!=0)) d <- max(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0)) 
              else {
                st <- 1 # if there is no valid dose below current dose
                earlystop = 1
              }
            }
            else {
              # when dec = "S", if there are valid doses above the current dose, escalate
              # if not, de-escalate, if neither, stop
              if (sum(which(dfut+dtox==0)<d)!=0) d <- max(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0))
              else {
                st <- 1
                earlystop = 1
              }
            }
          }
          # Below is the condition where the current dose is efficacious and not too toxic
          else {
            if (dec=="E") {
              # if there is a higher dose available, escalate
              if (d!=ndose && length(which(dtox[(d+1):ndose]+dfut[(d+1):ndose]==0))!=0) d <- min(which(dtox[(d+1):ndose]+dfut[(d+1):ndose]==0))+d 
              else d <- d
            }
            else if (dec=="D"){
              # if there is a lower dose available, de-escalate
              if (d!=1 && (length(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0)!=0))) d <- max(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0))  # else?
              else d <- d
            }
            else if (dec=="S"){
              d <- d
            }
          }
        }
      }
    }
    
    trial = simul 
    ntrial = nsimul
    
    
    if (earlystop==0){
      yT_c = x
      yE_c = y 
      elimi = dtox
      elimiE = ifelse(dfut == 2,1,0)
      
      pT_est<-(yT_c+0.05)/(n+0.1)
      pE_est<-(yE_c+0.05)/(n+0.1)
      pT_est<-pava(pT_est,n+0.1)+0.001*seq(1,ndose)
      
      pE_est<-peestimate(yE_c,n)			
      
      u<-u1*pE_est+(1-pT_est)*u2
      
      u[elimi==1]<--100
      u[elimiE==1]<--100
      u[n==0]<--100
      #u[pT_est>(targetT+0.1)]<--100
      
      d_mtd<-which.min(abs(pT_est-targetT))
      d_opT_est<-which.max(u[1:d_mtd])	
      dselect[trial]=d_opT_est
      
      if(d_opT_est==bd){bd.sel<-bd.sel+1/ntrial*100}
      if(abs(u.true[d_opT_est]-u.true[bd])<=(0.05*u.true[bd])&d_opT_est<=mtd){od.sel<-od.sel+1/ntrial*100}
      
	  if(pT.true[d_opt]>(targetT+0.1)){ov.sel<-ov.sel+1/ntrial*100}
	  
      dselect[trial]=d_opT_est
      
    } else {dselect[trial]<-99}	
    
    earlystop<-sum(dselect==99)/ntrial*100
    if(n[bd]<(npts/ndose)){poorall<-poorall+1/ntrial*100}
    overdose<-overdose+sum(n[pT.true>(targetT+0.1)])/ntrial/npts*100
    bd.pts<-bd.pts+n[bd]/ntrial/npts*100
    od.pts<-od.pts+sum(n[abs(u.true[1:mtd]-u.true[bd])<=(0.05*u.true[bd])])/ntrial/npts*100
    
    
    
  }
  
  
  results=list(bd.sel=bd.sel,od.sel=od.sel,bd.pts=bd.pts,od.pts=od.pts,
               earlystop=earlystop,ntox=0,neff=0,u.mean=0,
               overdose=overdose,poorall=poorall,incoherent=0,ov.sel=ov.sel)
  
  return(results)	
  
}


outputmat = NULL
for(utype in c(1,2)){
for(rtype in 1:1){
for(i in seq(8,30,by=2)){

oc =ji3plus3_sim(pT=0.3,  pE=0.4, eps1=0.05, eps2=0.05, 
                         psafe=0.95, pfutility=0.9, 
                         cohortsize=3,ncohort=i,startdose=1,
                         nsimul=10000,utilitytype = utype, randomtype = rtype)
						 

print(i)
outputmat=rbind(outputmat,c(i,utype,rtype,c(oc$bd.sel,oc$od.sel,oc$bd.pts,oc$od.pts,oc$earlystop,oc$overdose,oc$poorall)))
}}}


cname = c("ncohort","utype","rtype","bd.sel","od.sel","bd.pts","od.pts","earlystop","overdose","poorall")
colnames(outputmat)=cname
cname = c(cname,"design")
outputmat = cbind(outputmat,rep("ji3",nrow(outputmat)))
colnames(outputmat)=cname


write.csv(outputmat,paste0(iii,"/ji3_random.csv"),row.names=FALSE)

}

















