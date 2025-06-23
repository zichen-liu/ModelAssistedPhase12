
## ----------------------------------------------------------------------------
## Title: modified_boin12.R
## ----------------------------------------------------------------------------
## Authors: Haolun Shi, Ruitao Lin, Xiaolei Lin
## ----------------------------------------------------------------------------
## Description: Function to generate operating characteristics of BOIN12 
## ----------------------------------------------------------------------------
## Required Packages:  Iso, data.table
## ----------------------------------------------------------------------------
## Usage: get.oc.obd<-function(targetT=0.3,targetE=0.2,ncohort=10,cohortsize=3, startdose=1,
##                    cutoff.eli.T=0.95, cutoff.eli.E=0.95,u1,u2,ntrial=100,utilitytype=0,randomtype=1)  
##
##      targetT: target toxicity probability, i.e., the \phi_T in the manuscript.
##      targetE: lower limit for the efficacy rate, i.e., the \zeta_E in the manuscript.
##      ncohort: the number of cohorts.
##      cohortsize: the size of a cohort.
##      startdose: the starting dose level.
##      cutoff.eli.T: the posterior probability cutoff for toxicity dose elimination rule.
##      cutoff.eli.E: the posterior probability cutoff for futility dose elimination rule.
##      u1: the utility parameter w_11, in the scale of 0 to 100. 
##      u2: the utility parameter w_00, in the scale of 0 to 100. 
##      ntrial: the number of random trial replications.
##      utilitytype: a overriding argument for controlling the type of utility parameters. If set as 1, then (w_11,w_00) = (0.6,0.4). If set as 2, (w_11,w_00) = (1,0). Otherwise, the user-specified values for u1 and u2 are used. 
##      randomtype: not used.
##
## ----------------------------------------------------------------------------
## Value: 
##        bd.sel: OBD selection percentage
##        od.sel: favorable dose selection percentage
##        bd.pts: average percentage of patients at the OBD 
##        od.pts: average percentage of patients at the favorable doses 
##        earlystop: percentage of early stopped trials
##        overdose: overdose patients percentage 
##        poorall: poor allocation percentage
##        ov.sel: overdose selection percentage
##
## ----------------------------------------------------------------------------

library(Iso)
ndose <<- NDOSE
OBD <<- 0;

for (iii in 1:ndose){
OBD <<- iii 


simprob<-function(ndose,targetE,targetT,u1,u2 ,randomtype){
  targetE=TARGET.E 
  targetT=TARGET.T
  
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
  
  
#############################################################
#  Function to generate operating characteristics of BOIN12  #
#############################################################

get.oc.obd<-function(targetT=0.3,targetE=0.2,ncohort=10,cohortsize=3,
                    startdose=1,
                     cutoff.eli.T=0.95,
                      cutoff.eli.E=0.95,u1,u2,ntrial=100,utilitytype = 0,randomtype = 1){
  
 
  
if (utilitytype==1){
u1 = 60 
u2 = 40 
}
if (utilitytype==2){
u1 = 100 
u2 = 0 
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
  

n.earlystop=ncohort*cohortsize
  p.saf=0.6*targetT
  p.tox=1.4*targetT
  N1=6
  N2=9
  randomtype = 1
  ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
  pava <- function(x, wt = rep(1, length(x))) {
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
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }

  ## to make get.oc as self-contained function, we copied functions get.boundary() and select.mtd() here.
  get.boundary <- function(target, targetE, ncohort, cohortsize=3, p.saf=NA, p.tox=NA,  cutoff.eli=0.95,
                           cutoff.eli.E=0.90)
  {

    # if the user does not provide p.saf and p.tox, use the default values
    if(is.na(p.saf)) p.saf=0.6*target;
    if(is.na(p.tox)) p.tox=1.4*target;

    ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
    npts = ncohort*cohortsize;
    ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;elimE=NULL;
    for(n in (1:ncohort)*cohortsize)
    {
      error.min=3;
      for(m1 in 0:(n-1))
      {
        for(m2 in (m1+1):n)
        {

          error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
          error2 = 1-pbinom(m1, n, p.saf);
          error3 = pbinom(m2-1, n, p.tox);

          error=error1+error2+error3;
          if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
        }
      }
      ntrt = c(ntrt, n);
      b.e = c(b.e, cutoff1);
      b.d = c(b.d, cutoff2);

      elimineed=0; # indicating whether elimination is needed
      elimineedE=0
      if(n<3) { elim = c(elim, NA); elimE = c(elimE,NA)}  # require treating at least 3 patients before eliminating a dose
      else
      {
        for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
        {
          if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
        }
        if(elimineed==1) { elim = c(elim, ntox); }
        else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs

        for(neff in n:0){
          if(pbeta(targetE,neff+1,n-neff+1)>cutoff.eli.E){elimineedE=1; break;}
        }
        if(elimineedE==1){elimE=c(elimE,neff)} else {elimE=c(elimE,NA)}
      }
    }
    for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
    boundaries = rbind(ntrt, elim, b.d, b.e,elimE);
    rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=",
                             "Deescalate if # of DLT >=",  "Escalate if # of DLT <=", "Eliminate if # of Eff <=");
    colnames(boundaries) = rep("", ncohort);

    return(boundaries);
  }

  # function to compute Pr(Tox=1,Eff=1)
  f1<-function(x,bn.m1,bn.m2,rho){
    ff<-dnorm(x,bn.m1,1)*(1-pnorm(0,bn.m2+rho*(x-bn.m1),sqrt(1-rho^2)))
    return(ff)
  }



  poorall=0;
  incoherent=0;
  overdose=0;
  bd.sel=0;
  bd.pts=0;
  od.sel=0;
  ov.sel=0
  od.pts=0;
  ntox=0
  neff=0
  
  
  npts = ncohort*cohortsize;
  YT=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
  YE=matrix(rep(0, ndose*ntrial), ncol=ndose); # store efficacy outcome
  N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
  dselect = rep(0, ntrial); # store the selected dose level
  durationV = rep(0,ntrial)
  sel=rep(0,ndose);
  pts=rep(0,ndose);
  dlt=rep(0,ndose);
  eff=rep(0,ndose);
  poorall=0;
  temp=get.boundary(targetT, targetE, ncohort, cohortsize,cutoff.eli=cutoff.eli.T, cutoff.eli.E = cutoff.eli.E)
  b.e=temp[4,];   # escalation boundary
  b.d=temp[3,];   # deescalation boundary
  b.elim=temp[2,];  # elimination boundary
  b.elimE=temp[5,]
  u01=100
  u10=0
  u11 = u1 
  u00 = u2
  utility=c(u11,u10,u01,u00)
  # Assume independence between toxicity and efficacy
  targetP<-c(targetE*targetT,targetT*(1-targetE),(1-targetT)*targetE,(1-targetT)*(1-targetE))

  # Calculate the benchmark utility
  uu = sum(targetP*utility) # highest unacceptable utility is also the benchmark utility (i.e., desirable utility)

  # Calculate true utility
  #from marginal prob to joint prob
  get.joint.p<-function(pE,pT,r){
    p.ab<-function(pE,pT,a,b,r){
      pE^a*(1-pE)^(1-a)*pT^b*(1-pT)^(1-b)+(-1)^(a+b)*pE*(1-pE)*pT*(1-pT)*((exp(r)-1)/(exp(r)+1))
    }
    p=c(p.ab(pE,pT,a=1,b=0,r),p.ab(pE,pT,a=0,b=0,r),p.ab(pE,pT,a=1,b=1,r),p.ab(pE,pT,a=0,b=1,r))
    return(p)
  }

  p10<-p01<-p00<-p11<-rep(0,ndose)
if(FALSE){
  for(d in 1:ndose){
    p11[d]<-integrate(f1,bn.m1=qnorm(pT.true[d]),bn.m2=qnorm(pE.true[d]),rho=rho,lower=0,upper=Inf)$value
    p10[d]<-pT.true[d]-p11[d]
    p01[d]<-pE.true[d]-p11[d]
    p00[d]<-1-p11[d]-p10[d]-p01[d]
  }
}
  #u.true<-u11*p11+u00*p00+100*p01;
  
 set.seed(30); 
  for(trial in 1:ntrial){
    probs = simprob(ndose,targetE,targetT,u1,u2 ,randomtype)
    jj = probs$pE
	kk = probs$pT
		
	pE.true<-jj
	pT.true<-kk
	u.true<-(u1*pE.true+(1-pT.true)*u2);
    bd=probs$obd
	mtd = probs$mtd
    yT<-yE<-rep(0, ndose);    ## number of DLT/efficacy at each dose level
    y01<-y10<-y11<-y00<-rep(0,ndose); ## number of different outcomes at each dose level
    n<-rep(0, ndose);         ## number of patients treated at each dose level
    earlystop=0;              ## indicate whether the trial terminates early
    d=startdose;              ## starting dose level
    elimi = rep(0, ndose);    ## indicate whether doses are eliminated due to toxicity
    elimiE=  rep(0,ndose);    ## indicate whether doses are eliminated due to efficacy
    safe = 0
    posH<-rep(1-uu/100,ndose)
    duration=0
    for(i in 1:ncohort){

      

      T.time=0#obscohort$t.tox
      E.time=0#obscohort$t.tox
      inter.arrival=0#cumsum(rexp(cohortsize,rate=accrual.rate))
      t.all.seen=0#inter.arrival+pmax(T.time,E.time)
      duration=duration+max(t.all.seen)

     
      n[d]<-n[d]+cohortsize
	  
	  wT = sum(runif(cohortsize)<pT.true[d])
      yT[d] = yT[d] + wT;
      wE = sum(runif(cohortsize)<pE.true[d])
      yE[d] = yE[d] + wE;

      
      
      nc<-n[d]/cohortsize

      # determine whether current dose level is overly toxic
      if(!is.na(b.elim[nc]))
      {
        if(yT[d]>=b.elim[nc])
        {
          elimi[d:ndose]=1;
          if(d==1) {earlystop=1; break;}
        }
      }

      if(!is.na(b.elimE[nc]))
      {
        if(yE[d]<=b.elimE[nc])
        {
          elimi[d]=1;
        }
      }
      if(sum(elimi==1)==ndose) {earlystop=1; break;}

      u_curr<-(u1*yE[d]/n[d]+(1-yT[d]/n[d])*u2)/100*n[d] 	  
      posH[d] = 1-pbeta(uu/100,1+u_curr,n[d]-u_curr+1)
      posH <- posH*(1-elimi);
      if(n[d]>=N1){safe=1} else{safe=0}
      if(n[d]>=n.earlystop){break}

      if (yT[d]>=b.d[nc] && d!=1) {
        if(sum(elimi[1:(d-1)]==0)>0){d_opt=max(which(elimi[1:(d-1)]==0))} else {
          if(elimi[d]==1){earlystop=1;break} else{d_opt=d}}
      } else if (yT[d]>=b.d[nc] && d==1) {if(elimi[d]==0){d_opt=d} else{earlystop=1;break}
      } else{
        admi_set=d;
        if(d>1){
          if(sum(elimi[1:(d-1)]==0)>0){admi_set<-c(admi_set,max(which(elimi[1:(d-1)]==0)))}
        }
        if(d<ndose){
          if(safe==0){
            if(sum(elimi[(d+1):ndose]==0)>0){admi_set<-c(admi_set,d+min(which(elimi[(d+1):ndose]==0)))}
          } else {
            if(yT[d]<=b.e[nc] & sum(elimi[(d+1):ndose]==0)>0){admi_set<-c(admi_set,d+min(which(elimi[(d+1):ndose]==0)))}
          }
        }
        #if(length(admi_set)>1 & eff_cut>0.5 & n[d]>=9){admi_set<-admi_set[-1]}
        temp.posH<-posH[admi_set]+runif(length(admi_set))*(10^-15)
        d_opt=admi_set[which.max(temp.posH)]
      }

      if (elimi[d_opt]==1) {earlystop=1; break}
      if (sum(elimi)==ndose) {earlystop=1; break}

      if (d<ndose){
        if(sum(elimi[(d+1):ndose]==0)>0){
          d_temp=d+min(which(elimi[(d+1):ndose]==0))
          if(n[d]>=N2 & n[min(d_temp,ndose)]==0 & yT[d]<b.d[n[d]/cohortsize]){d_opt<-d_temp}
        }
      }

      d<-d_opt
    }
    YT[trial,]=yT;
    YE[trial,]=yE;
    N[trial,]=n;
    durationV[trial]=duration
    if (earlystop==0){
	
	
	  
	  pT<-(yT+0.05)/(n+0.1)
      pE<-(yE+0.05)/(n+0.1)
      pT<-pava(pT,n+0.1)+0.001*seq(1,ndose)
      pE<-peestimate(yE,n)
      u<-u1*pE+(1-pT)*u2
      u[elimi==1]<--100
      u[elimiE==1]<--100
      u[n==0]<--100
      #u[pT>(targetT+0.1)]<--100
      d_mtd<-which.min(abs(pT-targetT))
      d_opt<-which.max(u[1:d_mtd])	
      dselect[trial]=d_opt
	  
	  if(d_opt==bd){bd.sel<-bd.sel+1/ntrial*100}
      if(abs(u.true[d_opt]-u.true[bd])<=(0.05*u.true[bd]) & d_opt<=mtd){od.sel<-od.sel+1/ntrial*100}
	  if(pT.true[d_opt]>(targetT+0.1)){ov.sel<-ov.sel+1/ntrial*100}
      sel[dselect[trial]]=sel[dselect[trial]]+1/ntrial*100
	  
    } else {dselect[trial]<-99}
    
	earlystop<-sum(dselect==99)/ntrial*100
    if(n[bd]<(npts/ndose)){poorall<-poorall+1/ntrial*100}
    overdose<-overdose+sum(n[pT.true>(targetT+0.1)])/ntrial/npts*100
    bd.pts<-bd.pts+n[bd]/ntrial/npts*100
    od.pts<-od.pts+sum(n[abs(u.true[1:mtd]-u.true[bd])<=(0.05*u.true[bd])])/ntrial/npts*100
	
	pts<-pts+n/ntrial
    dlt<-dlt+yT/ntrial
    eff<-eff+yE/ntrial
	
	ntox<-ntox+sum(yT)/ntrial
    neff<-neff+sum(yE)/ntrial
	
  }
  sel<-round(sel,1)
  pts<-round(pts,1)
  dlt<-round(dlt,1)
  u.true<-round(u.true,1)
  earlystop<-sum(dselect==99)/ntrial*100
  
  results=list(bd.sel=bd.sel,od.sel=od.sel,bd.pts=bd.pts,od.pts=od.pts,
               earlystop=earlystop,ntox=ntox,neff=neff,u.mean=0,
               overdose=overdose,poorall=poorall,incoherent=0,ov.sel=ov.sel)
			   
	return(results)

}





outputmat = NULL
for(utype in c(1,2)){
for(rtype in c(1)){
for(i in SSIZERANGE){



oc = get.oc.obd(targetT=TARGET.T,targetE=TARGET.E,ncohort=i,cohortsize=3,
                    startdose=STARTD,
                     cutoff.eli.T=0.95,
                     cutoff.eli.E=0.9,u1,u2,ntrial=10000,utilitytype = utype,randomtype = rtype)
					 
print(i)
outputmat=rbind(outputmat,c(i,utype,rtype,c(oc$bd.sel,oc$od.sel,oc$bd.pts,oc$od.pts,oc$earlystop,oc$overdose,oc$poorall,oc$ov.sel)))
}}}



cname = c("ncohort","utype","rtype","bd.sel","od.sel","bd.pts","od.pts","earlystop","overdose","poorall","ov.sel")
colnames(outputmat)=cname
cname = c(cname,"design")
outputmat = cbind(outputmat,rep("modified_boin12",nrow(outputmat)))
colnames(outputmat)=cname


write.csv(outputmat,paste0(iii,"/modified_boin12_random.csv"),row.names=FALSE)

}
