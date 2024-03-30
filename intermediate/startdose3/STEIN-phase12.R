setwd('/tmp/temp')


library(Iso)



OBD <<- 0;

for (iii in 1:5){
OBD <<- iii 


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
  
get.oc <- function(targetT,  targetE, psi1, psi2, ncohort, cohortsize, startdose=1, cutoff.eli=0.95, ntrial=10, utilitytype = 1, randomtype = 1)
{
target = targetT

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

	## to make get.oc as self-contained function, we copied functions get.boundary() and select.mtd() here.
	get.boundary <- function(target, ncohort, cohortsize=3, p.saf=NA, p.tox=NA,  cutoff.eli=0.95)
	{
 
# if the user does not provide p.saf and p.tox, use the default values
		if(is.na(p.saf)) p.saf=0.75*target;
		if(is.na(p.tox)) p.tox=1.25*target;

### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
		npts = ncohort*cohortsize;
		ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;
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
			if(n<3) { elim = c(elim, NA); }  # require treating at least 3 patients before eliminating a dose
			else
			{
				for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
				{
					if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
				}
				if(elimineed==1) { elim = c(elim, ntox); }
				else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
			}
		}
		for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
		boundaries = rbind(ntrt, elim, b.d, b.e);
		rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=", 
								 "Deescalate if # of DLT >=",  "Escalate if # of DLT <=");
		colnames(boundaries) = rep("", ncohort);
		
		return(boundaries);
	}
	


ndose = 5
npts = ncohort*cohortsize;
YT=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
YE=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
dselect = rep(0, ntrial); # store the selected dose level
sel=rep(0,ndose);
pts=rep(0,ndose);
dlt=rep(0,ndose);
eff=rep(0,ndose);
ntox=0;
neff=0;
acr=0;
exc=0;
temp=get.boundary(targetT, ncohort, cohortsize,cutoff.eli=0.95) 	
b.e=temp[4,];   # escalation boundary
b.d=temp[3,];   # deescalation boundary
b.elim=temp[2,];  # elimination boundary
psi<-log((1-psi1)/(1-psi2))/log(psi2*(1-psi1)/psi1/(1-psi2))

bd.sel=0;
bd.pts=0;
od.sel=0;
od.pts=0;
ntox=0;
neff=0;
poorall=0;
incoherent=0;
overdose=0;


#print(cutoff.e)
set.seed(30);
	################## simulate trials ###################
	for(trial in 1:ntrial)
	{

		yT<-yE<-rep(0, ndose);    ## the number of DLT at each dose level
		n<-rep(0, ndose);    ## the number of patients treated at each dose level
		earlystop=0;         ## indiate whether the trial terminates early
		d=startdose;         ## starting dose level
		elimi = rep(0, ndose);  ## indicate whether doses are eliminated
		elimiE=  rep(0,ndose);
		
		probs = simprob(ndose,targetE,targetT,u1,u2 ,randomtype)
		jj = probs$pE
		kk = probs$pT
		pE.true<-jj
		pT.true<-kk
		u.true<-(u1*pE.true+(1-pT.true)*u2);
		bd  =probs$obd
		mtd =probs$mtd
		for(i in 1:ncohort)  
		{  			
		
		    
			### generate toxicity outcome
			
			wT = sum(runif(cohortsize)<pT.true[d])
			yT[d] = yT[d] + wT;
			wE = sum(runif(cohortsize)<pE.true[d])
			yE[d] = yE[d] + wE;
			n[d] = n[d] + cohortsize;
			nc = n[d]/cohortsize;

			if(!is.na(b.elim[nc]))
			{
				if(yT[d]>=b.elim[nc]) 
				{      
					elimi[d:ndose]=1;
					if(d==1) {earlystop=1; break;} 
				}
				
			}
			if(n[d]>=3 && pbeta(targetE,yE[d]+1,n[d]-yE[d]+1)>0.9) {elimiE[d]=1;}
			
 

			if (yT[d]>=b.d[nc] && d!=1) {d_opt=max(which(elimi[1:(d-1)]==0))
			} else if (yT[d]>=b.d[nc] && d==1) {d_opt=d
			} else{
			       if (yE[d]/n[d]>psi) {d_opt=d
			} else if (yT[d]> b.e[nc]) {
			posH<-c(pbeta(psi,yE[max(1,d-1)]+1,n[max(1,d-1)]-yE[max(1,d-1)]+1),
			pbeta(psi,yE[d]+1,n[d]-yE[d]+1))-c(0,0.001)
            d_opt<-max(d-2+which.min(posH),1);
			} else if (yT[d]<=b.e[nc] && d==1) {
			       if (elimi[d+1]==0) {
			posH<-c(pbeta(psi,yE[d]+1,n[d]-yE[d]+1),
			(n[d+1]>0)*pbeta(psi,yE[d+1]+1,n[d+1]-yE[d+1]+1))-c(0,0.001)	
			d_opt=d-1+which.min(posH)} else {d_opt=d}
			} else if (yT[d]<=b.e[nc] && d==ndose) {
			posH<-c(pbeta(psi,yE[d-1]+1,n[d-1]-yE[d-1]+1),
			pbeta(psi,yE[d]+1,n[d]-yE[d]+1))-c(0,0.001)
			d_opt=d-2+which.min(posH)
		    } else if (yT[d]<=b.e[nc] && d!=1 && d!=ndose) {
		    	       if (elimi[d+1]==0){
		    	posH<-c(pbeta(psi,yE[d-1]+1,n[d-1]-yE[d-1]+1),pbeta(psi,yE[d]+1,n[d]-yE[d]+1),
				(n[d+1]>0)*pbeta(psi,yE[d+1]+1,n[d+1]-yE[d+1]+1))-c(0,0.001,0.002)       	
				d_opt=d-2+which.min(posH)		    	       	
		    	       } else {
			posH<-c(pbeta(psi,yE[d-1]+1,n[d-1]-yE[d-1]+1),
			pbeta(psi,yE[d]+1,n[d]-yE[d]+1))-c(0,0.001)
			d_opt=d-2+which.min(posH)		    	       	
		    	       	
		    	       }

			}
            }	
            if (elimiE[d_opt]==1) {earlystop=1; break} 
			d<-d_opt


		}

		YT[trial,]=yT;
		YE[trial,]=yE;
		N[trial,]=n;
        ntox=ntox+sum(yT)/ntrial
		neff=neff+sum(yE)/ntrial
		if (earlystop==0){
		pT<-(yT+0.05)/(n+0.1)
		pE<-(yE+0.05)/(n+0.1)
		pT<-pava(pT,n+0.1)+0.001*seq(1,ndose)

        pE<-peestimate(yE,n)
		
		
		if(FALSE){  ##Old way of finding the OBD
		pT[elimi==1]<-20
		pT[n==0]<-20
		pE[n==0]<-0
		pE[elimi==1]<-0
		u<-pE-0.33*pT-1.09*(pT>targetT) 	
		d_opt<-which.max(u)	
		}
		
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
	  
		dselect[trial]=d_opt
		sel[dselect[trial]]=sel[dselect[trial]]+1/ntrial*100
		#if (d_opt==2) {print(yT);print(yE);print(n)}
		} else {dselect[trial]<-99}	
		
		pts<-pts+n/ntrial
		dlt<-dlt+yT/ntrial
		eff<-eff+yE/ntrial	
		
		earlystop<-sum(dselect==99)/ntrial*100
		if(n[bd]<(npts/ndose)){poorall<-poorall+1/ntrial*100}
		overdose<-overdose+sum(n[pT.true>(targetT+0.1)])/ntrial/npts*100
		bd.pts<-bd.pts+n[bd]/ntrial/npts*100
		od.pts<-od.pts+sum(n[abs(u.true[1:mtd]-u.true[bd])<=(0.05*u.true[bd])])/ntrial/npts*100
				
	}	
	pts<-round(pts,1)
	dlt<-round(dlt,1)
	
	# results:
	# pT.true: true toxicity rate
	# pE.true: true efficacy rate
	# sel: selection percentage of each dose
	# pts: number of patients treated at each dose
	# dlt: number of DLTs observed at each dose
	# eff: number of efficacy outcomes observed at each dose
	# ntox: total number of DLTs
	# neff: total number of Effs
	results=list(bd.sel=bd.sel,od.sel=od.sel,bd.pts=bd.pts,od.pts=od.pts,
               earlystop=earlystop,ntox=ntox,neff=neff,u.mean=0,
               overdose=overdose,poorall=poorall,incoherent=0)
			   
	#results=list(pT.true=pT.true,pE.true=pE.true,sel=sel,pts=pts,dlt=dlt,eff=eff,ntox=ntox,neff=neff)
	return(results)	
}



cohortsize = 3 # cohort size
ncohort = 10 # number of cohorts
psi1<-0.2 # lowest acceptable efficacy rate
psi2<-0.6 # very promising efficacy rate that leads to dose retainment





outputmat = NULL
for(utype in c(1,2)){
for(rtype in 1:1){
for(i in c(12,18)){

oc = get.oc(targetT=0.3,  targetE=0.2, psi1, psi2, ncohort=i, cohortsize=3, startdose=3, cutoff.eli=0.95, ntrial=10000, utilitytype = utype, randomtype = rtype)


print(i)
outputmat=rbind(outputmat,c(i,utype,rtype,c(oc$bd.sel,oc$od.sel,oc$bd.pts,oc$od.pts,oc$earlystop,oc$overdose,oc$poorall)))
}}}


cname = c("ncohort","utype","rtype","bd.sel","od.sel","bd.pts","od.pts","earlystop","overdose","poorall")
colnames(outputmat)=cname
cname = c(cname,"design")
outputmat = cbind(outputmat,rep("stein",nrow(outputmat)))
colnames(outputmat)=cname


write.csv(outputmat,paste0(iii,"/stein_random.csv"),row.names=FALSE)


}