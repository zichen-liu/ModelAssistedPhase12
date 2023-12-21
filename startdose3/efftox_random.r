setwd('/tmp/temp')


OBD <<- 0;

for (iii in 1:5){
OBD <<- iii 

ncohort = 10

ntrial=2;
library(trialr)

efftox<-function(ncohort,ntrial,utilitytype = 1, randomtype = 1){


if (utilitytype==1){
u1 = 60 
u2 = 40 
}
if (utilitytype==2){
u1 = 100 
u2 = 0 
}


targetT = 0.3
targetE = 0.2
cohortsize = 3

npts = cohortsize*ncohort
uu = u1*0.7+(1-0.1)*u2
cutoff.eliT=0.95
cutoff.eliE=0.9
ndose = 5
res=NULL;
dselect = rep(0, ntrial); # store the selected dose level
bd.sel=0;
bd.pts=0;
od.sel=0;
od.pts=0;
ntox=0;
neff=0;
poorall=0;
incoherent=0;
overdose=0;
u.mean=0;

if(utilitytype==1){
eff0 = 0.20
  tox1 = 0.533
  eff_star = 0.5
  tox_star = 0.20
}

if(utilitytype==2){
eff0 = 0.80
  tox1 = 0.99
  eff_star = 0.90
  tox_star = 0.50
}

  
pp <- efftox_solve_p(eff0 = eff0, tox1 = tox1, eff_star = eff_star, tox_star =tox_star)

dat <- list(
  num_doses = 5,
  real_doses = c(1, 2, 3, 4, 5),
  efficacy_hurdle = 0.2,
  toxicity_hurdle = 0.3,
  p_e = 0.05,
  p_t = 0.1,
  p = pp,
 eff0 = eff0, tox1 = tox1, eff_star = eff_star, tox_star =tox_star,  
  alpha_mean = -2.823, alpha_sd = 2.7099,
  beta_mean = 3.9364, beta_sd = 2.7043,
  gamma_mean = -2.8240, gamma_sd = 2.7108,
  zeta_mean = 3.9374, zeta_sd = 2.7032,
  eta_mean = 0, eta_sd = 0.2,
  psi_mean = 0, psi_sd = 1,
  
  doses = c(),
  tox   = c(),
  eff   = c(),
  num_patients = 0
)


set.seed(30);
for(trial in 1:ntrial){
  
  
  
	probs = simprob(ndose,targetE,targetT,u1,u2 ,randomtype)
	
	jj = probs$pE
kk = probs$pT
	
    pE.true<-jj
    pT.true<-kk
    u.true<-(u1*pE.true+(1-pT.true)*u2);
    bd  =probs$obd
    mtd  =probs$mtd
	
  
  temp = efftox_simulate(dat, num_sims = 1, first_dose = 3, 
                         true_eff = pE.true,
                         true_tox = pT.true,
                         cohort_sizes = rep(cohortsize, ncohort),chains=1,iter = 500,show_messages = FALSE,open_progress=FALSE,refresh = 0)

  if(is.na(temp$recommended_dose)){dselect[trial]=99} else {
    dselect[trial] = temp$recommended_dose
  }
  if(dselect[trial]<99){
    if(dselect[trial]==bd){bd.sel<-bd.sel+1/ntrial*100}
	d_opt = dselect[trial]		   
	
	if(abs(u.true[d_opt]-u.true[bd])<=(0.05*u.true[bd]) & d_opt<=mtd){od.sel<-od.sel+1/ntrial*100}
  }
  earlystop<-sum(dselect==99)/ntrial*100
  n<-rep(0,ndose)
  yE<-rep(0,ndose)
  yT<-rep(0,ndose)
  for(j in 1:ndose){
    n[j]<-sum(as.numeric(unlist(temp$doses_given))==j)
    yE[j]<-sum(as.numeric(unlist(temp$efficacies))[as.numeric(unlist(temp$doses_given))==j])
    yT[j]<-sum(as.numeric(unlist(temp$toxicities))[as.numeric(unlist(temp$doses_given))==j])
  }
  bd.pts<-bd.pts+n[bd]/ntrial/npts*100
  od.pts<-od.pts+sum(n[abs(u.true[1:mtd]-u.true[bd])<=(0.05*u.true[bd])])/ntrial/npts*100  
  if(n[bd]<(npts/ndose)){poorall<-poorall+1/ntrial*100}
  overdose<-overdose+sum(n[pT.true>(targetT+0.1)])/ntrial/npts*100
}



  results=list(bd.sel=bd.sel,od.sel=od.sel,bd.pts=bd.pts,od.pts=od.pts,
               earlystop=earlystop,ntox=ntox,neff=neff,u.mean=u.mean,
               overdose=overdose,poorall=poorall,incoherent=0)
  return(results)	
}









xmat = NULL
for(utype in c(1,2)){
for(rtype in c(1)){
for(i in c(12,18)){
xmat = rbind(xmat,c(i,utype,rtype))
}}}



library(doMC)
library(doRNG)
registerDoMC(64)


outputmat = NULL

#ooo = foreach(trial=1:ntrial,.combine=rbind,.errorhandling='pass')%dopar% {
ooo = foreach(kk=1:nrow(xmat),.combine=rbind)%dorng% {
	
i = xmat[kk,1]
utype = xmat[kk,2]
rtype=xmat[kk,3]

oc= efftox(ncohort=i,ntrial=500,utilitytype = utype, randomtype = rtype)

print(i)
return(c(i,utype,rtype,c(oc$bd.sel,oc$od.sel,oc$bd.pts,oc$od.pts,oc$earlystop,oc$overdose,oc$poorall)))
}

outputmat = ooo 
cname = c("ncohort","utype","rtype","bd.sel","od.sel","bd.pts","od.pts","earlystop","overdose","poorall")
colnames(outputmat)=cname
cname = c(cname,"design")
outputmat = cbind(outputmat,rep("efftox",nrow(outputmat)))
colnames(outputmat)=cname


write.csv(outputmat,paste0(iii,"/efftox_random.csv"),row.names=FALSE)


}