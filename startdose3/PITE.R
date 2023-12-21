setwd('/tmp/temp')




OBD <<- 0;

for (iii in 1:5){
OBD <<- iii 


###############################################################
##########   simulate nsimul trials   #########################
PITE_sim <- function(pT, pE, psafe=0.95, pfutility=0.95, eps1=0.05, eps2=0.05, startdose=1, cohortsize=3,ncohort=10, nsimul=1000,utilitytype = 1, randomtype = 1){


if (utilitytype==1){
u1 = 60 
u2 = 40 
}
if (utilitytype==2){
u1 = 100 
u2 = 0 
}


#########################################
# define parameter values
#########################################
# Calculate JUPM
upmgen <- function(tbound,ebound,parabeta,n,ntox,neff){
  upm <- matrix(0,length(tbound)-1,length(ebound)-1)
  for (i in 1:(length(tbound)-1)){
    for (j in 1:(length(ebound)-1)){
      upm[i,j]=( pbeta(tbound[i+1],parabeta[1]+ntox,parabeta[2]+n-ntox) - 
                   pbeta(tbound[i],parabeta[1]+ntox,parabeta[2]+n-ntox) ) * 
        ( pbeta(ebound[j+1],parabeta[1]+neff,parabeta[2]+n-neff) - 
            pbeta(ebound[j],parabeta[1]+neff,parabeta[2]+n-neff) ) / 
        ((tbound[i+1]-tbound[i])*(ebound[j+1]-ebound[j]))
    }
  }
  return(upm)
}

maptoDEC <- function(index){
  if (index %in% c("UU", "UL")){
    DEC <- "D"
  }else if(index %in% c("EL", "EU", "LU")){
    DEC <- "S"
  }else{DEC <- "E"}
  return(DEC)
}

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


  
  
#pT <- 0.3
qE <- 0.2


parabeta <- c(1,1) #Beta distribution parameters



#########################################
# create upper and lower interval for toxicity
up.int.bound.tox <- seq(pT+eps2, 1, min(eps1+eps2,1-pT-eps2))
mid.int.bound.tox <- c(pT-eps1, pT+eps2)
low.int.bound.tox <- sort(-seq(-pT+eps1, 0, min(eps1+eps2,pT-eps1)))
tbound <- c(0, low.int.bound.tox,up.int.bound.tox, 1)

up.int.bound.eff <- seq(pE, 1, 0.2)
low.int.bound.eff <- sort(-seq(-pE, 0, 0.2))
ebound <- c(low.int.bound.eff[-length(low.int.bound.eff)], up.int.bound.eff)




#########################################
# define the true decision table
maptoINT <- matrix(c(rep(c(rep("LL", length(low.int.bound.tox)), "EL",
                           rep("UL", length(up.int.bound.tox))), 2),
                     rep(c(rep("LU", length(low.int.bound.tox)),  "EU",
                           rep("UU", length(up.int.bound.tox))), 3)), 
                   ncol = 5, byrow = FALSE)

  


  csize= cohortsize
  sampsize = csize * ncohort 
  
  targetT = pT
  targetE = pE
  ntrial = nsimul
  npts = ncohort*cohortsize;
	  
	bd.sel=0;
	bd.pts=0;
	od.sel=0;
	od.pts=0;
	ntox=0;
	neff=0;
	poorall=0;
	incoherent=0;
	overdose=0;
  ndose <- 5

  dselect = rep(0, ntrial); # store the selected dose level
  set.seed(30);
  # run 1000 simulation
  for (simul in 1:nsimul){ 
  

	
	probs = simprob(ndose,targetE,targetT,u1,u2 ,randomtype)
	jj = probs$pE
	kk = probs$pT
	pE.true<-jj
	pT.true<-kk
	u.true<-(u1*pE.true+(1-pT.true)*u2);
	bd  =probs$obd
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
      #cat(simul, " ", file = file1, append = T, sep = ",")
      #cat(d, " ", file = file1, append = T, sep = ",")
	  
      xx <-sum(runif(cohortsize)<toxtrue[d]) # generate tox outcome at current dose
	  
	  #cat(xx, " ", file = file1, append = T, sep = ",")
      
	  yy <- sum(runif(cohortsize)<efftrue[d])
	   # generate eff outcome at current dose
      
	  #cat(yy, " ", file = file1, append = T, sep = ",")
      
      x[d] <- x[d] + xx # total tox outcome
      y[d] <- y[d] + yy # total eff outcome
      n[d] <- n[d] + csize # total sample size used
      currentsize <- currentsize + csize
      
      # safety rule	(monotonic tox)	
      if (pbeta(pT, parabeta[1] + x[d], parabeta[2] + n[d] - x[d]) < 1 - psafe) {
        dtox[d:ndose] <- 1
      }
      
      # futility rule (mark only current dose)
      if (pbeta(qE, parabeta[1] + y[d], parabeta[2] + n[d] - y[d]) > pfutility) {
        dfut[d] <- 2
      }
      
      upm <- upmgen(tbound, ebound, parabeta, n[d], x[d], y[d]) # compute JUPM
      index <- maptoINT[which(upm == max(upm),arr.ind=T)] # map to interval
      dec <- maptoDEC(index) # get decision
      #print(paste("current dose", d))
      #print(paste("current decision", dec))
      # if the sampsize is reached, stop
      if (currentsize >= sampsize){
        st <- 1
		
        #print("stop due to sample size")
      }
      else {
        # if all doses are either too toxic or of no efficacy (couldn't find a dose for the next corhort)
        if (length(which(dfut+dtox==0))==0) {
          st <- 1 
          #print("stop since all doses are either too toxic or of no efficacy")
		  earlystop = 1
        }
        else {
          # if current dose is too toxic, de-escalate
          if (dtox[d]==1) d <- max(which(dfut[1:(d-1)]==0))
          # if current dose is of no efficacy
          else if (dfut[d]==2) { 
            if (dec=="E") { 
              # if there is a valid dose above current dose, escalate by minimum dose size
              if (d!=ndose && length(which(dtox[(d+1):ndose]+dfut[(d+1):ndose]==0))!=0) d <- min(which(dtox[(d+1):ndose]+dfut[(d+1):ndose]==0))+d 
              # if there is no dose above current dose, de-escalate to the next available dose
              else d <- max(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0)) 
            }
            else if (dec=="D") {
              # if there is a valid dose below current dose, de-escalate by minimum dose size
              if (d!=1 && length(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0)!=0)) d <- max(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0)) 
              else {
                st <- 1 # if there is no valid dose below current dose
                #print("current dose not toxic but no efficacy, decision D")
				earlystop = 1
              }
            }
            else {
              # when dec = "S", if there are valid doses above the current dose, escalate
              # if not, de-escalate, if neither, stop
              if (sum(which(dfut+dtox==0)<d)!=0) d <- max(which(dtox[1:(d-1)]+dfut[1:(d-1)]==0))
              else {
                st <- 1
                #print("current dose not toxic but no efficacy, decision S")
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
      #print(paste("next dose", d))
      #cat("\n",file=file1,append=T)
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
        if(abs(u.true[d_opT_est]-u.true[bd])<=(0.05*u.true[bd]) & d_opT_est<=mtd){od.sel<-od.sel+1/ntrial*100}
	  
		dselect[trial]=d_opT_est
		#sel[dselect[trial]]=sel[dselect[trial]]+1/ntrial*100
		#if (d_opT_est==2) {print(yT_c);print(yE_c);print(n)}
	} else {dselect[trial]<-99}	
		
	
	earlystop<-sum(dselect==99)/ntrial*100
	if(n[bd]<(npts/ndose)){poorall<-poorall+1/ntrial*100}
	overdose<-overdose+sum(n[pT.true>(targetT+0.1)])/ntrial/npts*100
	bd.pts<-bd.pts+n[bd]/ntrial/npts*100
	od.pts<-od.pts+sum(n[abs(u.true[1:mtd]-u.true[bd])<=(0.05*u.true[bd])])/ntrial/npts*100
			
  }
  results=list(bd.sel=bd.sel,od.sel=od.sel,bd.pts=bd.pts,od.pts=od.pts,
               earlystop=earlystop,ntox=0,neff=0,u.mean=0,
               overdose=overdose,poorall=poorall,incoherent=0)
			   	
	return(results)	
				
  
  
}


###############################################################



outputmat = NULL
for(utype in c(1,2)){
for(rtype in 1:1){
for(i in c(12,18)){



oc = PITE_sim(pT=0.3, pE=0.4, psafe=0.95, pfutility=0.9, eps1=0.05, eps2=0.05,  startdose=3, cohortsize=3,ncohort=i, nsimul=10000,utilitytype = utype, randomtype = rtype)

print(i)
outputmat=rbind(outputmat,c(i,utype,rtype,c(oc$bd.sel,oc$od.sel,oc$bd.pts,oc$od.pts,oc$earlystop,oc$overdose,oc$poorall)))
}}}


cname = c("ncohort","utype","rtype","bd.sel","od.sel","bd.pts","od.pts","earlystop","overdose","poorall")
colnames(outputmat)=cname
cname = c(cname,"design")
outputmat = cbind(outputmat,rep("printe",nrow(outputmat)))
colnames(outputmat)=cname

write.csv(outputmat,paste0(iii,"/printe_random.csv"),row.names=FALSE)


}














