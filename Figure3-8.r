library(data.table)


setwd(paste0(PATH,'/intermediate/NDOSE5/'))

par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)


for(UTYPE in c(1 ,2)){

RTYPE = 1
getvalue<-function(path){
for(i in 1:5){
utpi= fread(paste0(i,path))
KK = ncol(utpi)-1
if (i==1){
utpix = utpi[,1:KK]
}else{
utpix = utpi[,1:KK] + utpix
}
}
utpix = utpix / 5
return(cbind(utpix,utpi$design))
}

utpi= getvalue("/utpi_random.csv")
boin=getvalue("/BOINET_random.csv")
efft=getvalue("/efftox_random.csv")
tepi=getvalue("/tepi_random.csv")
boin12 = getvalue("/boin12_random.csv")
modified_boin12 = getvalue("/modified_boin12_random.csv")
ji3 = getvalue("/ji3_random.csv")
printe = getvalue("/printe_random.csv")
stein = getvalue("/stein_random.csv")

tblall = NULL
tblall = rbind(tblall, utpi)
tblall = rbind(tblall, boin)
tblall = rbind(tblall, efft)
tblall = rbind(tblall, tepi)
tblall = rbind(tblall, boin12)
tblall = rbind(tblall, modified_boin12)
tblall = rbind(tblall, ji3)
tblall = rbind(tblall, printe)
tblall = rbind(tblall, stein)




titlep = paste(ifelse(UTYPE==1,"Utility OBD","Efficacy OBD"))

index = which(utpi$rtype==RTYPE & utpi$utype==UTYPE)


utpi= utpi[index,]
boin= boin[index,]
efft= efft[which(efft$rtype==RTYPE & efft$utype==UTYPE),]
tepi= tepi[index,]
boin12= boin12[index,]
modified_boin12= modified_boin12[index,]
ji3= ji3[index,]
printe= printe[index,]
stein= stein[index,]




sm = smooth.spline(efft$ncohort, efft$bd.sel, spar = 0.6)
efft$bd.sel =predict(sm,efft$ncohort)$y	
 
sm = smooth.spline(efft$ncohort, efft$od.sel, spar = 0.6)
efft$od.sel =predict(sm,efft$ncohort)$y	

sm = smooth.spline(efft$ncohort, efft$bd.pts, spar = 0.6)
efft$bd.pts =predict(sm,efft$ncohort)$y	 

sm = smooth.spline(efft$ncohort, efft$od.pts, spar = 0.6)
efft$od.pts =predict(sm,efft$ncohort)$y	 

sm = smooth.spline(efft$ncohort, efft$poorall, spar = 0.3)
efft$poorall =predict(sm,efft$ncohort)$y	 

sm = smooth.spline(efft$ncohort, efft$overdose, spar = 0.6)
efft$overdose =predict(sm,efft$ncohort)$y	 

sm = smooth.spline(efft$ncohort, efft$ov.sel, spar = 0.5)
efft$ov.sel =predict(sm,efft$ncohort)$y	 



ltys=seq(1,9)
liness=seq(1,9)
liness=rep(1,9)





## Figure 3
cs = c("Black","Red","Green","Blue",5,6,"Pink","Orange","Purple")

dname=c("uTPI", "PRINTE","BOIN-ET","TEPI","Joint3+3","BOIN12","STEIN","EffTox", "M-BOIN12")

# setEPS()
# postscript(paste0(PATH,"/results/Figure3_",UTYPE,".eps"))
png(paste0(PATH,"/results/Figure3_",UTYPE,".png"), width = 7, height = 7, units = "in", res = 300)

xx = c(utpi$bd.sel,
printe$bd.sel,
boin$bd.sel,
tepi$bd.sel,
efft$bd.sel,
boin12$bd.sel,
modified_boin12$bd.sel,
stein$bd.sel,
ji3$bd.sel)



par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)
plot(utpi$ncohort*3,utpi$bd.sel,ylab="Best dose selection %",xlab="Number of patients",type="l",pch=ltys[1],col=cs[1],xlim=c(24,90),ylim=c(min(xx)*0.95,max(xx)*1.05),lwd=1.5,lty=liness[1],main=titlep)
points(utpi$ncohort*3,utpi$bd.sel,type="p",pch=ltys[1],col=cs[1],lwd=1.5)

lines(printe$ncohort*3,printe$bd.sel,type="l",pch=ltys[2],col=cs[2],lwd=1.5,lty=liness[2])
points(printe$ncohort*3,printe$bd.sel,type="p",pch=ltys[2],col=cs[2],lwd=1.5)

lines(boin$ncohort*3,boin$bd.sel,type="l",pch=ltys[3],col=cs[3],lwd=1.5,lty=liness[3])
points(boin$ncohort*3,boin$bd.sel,type="p",pch=ltys[3],col=cs[3],lwd=1.5)

lines(tepi$ncohort*3,tepi$bd.sel,type="l",pch=ltys[4],col=cs[4],lwd=1.5,lty=liness[4])
points(tepi$ncohort*3,tepi$bd.sel,type="p",pch=ltys[4],col=cs[4],lwd=1.5)

lines(ji3$ncohort*3,ji3$bd.sel,type="l",pch=ltys[5],col=cs[5],lwd=1.5,lty=liness[5])
points(ji3$ncohort*3,ji3$bd.sel,type="p",pch=ltys[5],col=cs[5],lwd=1.5)

lines(boin12$ncohort*3,boin12$bd.sel,type="l",pch=ltys[6],col=cs[6],lwd=1.5,lty=liness[6])
points(boin12$ncohort*3,boin12$bd.sel,type="p",pch=ltys[6],col=cs[6],lwd=1.5)

lines(stein$ncohort*3,stein$bd.sel,type="l",pch=ltys[7],col=cs[7],lwd=1.5,lty=liness[7])
points(stein$ncohort*3,stein$bd.sel,type="p",pch=ltys[7],col=cs[7],lwd=1.5)

lines(efft$ncohort*3,efft$bd.sel,type="l",pch=ltys[8],col=cs[8],lwd=1.5,lty=liness[8])
points(efft$ncohort*3,efft$bd.sel,type="p",pch=ltys[8],col=cs[8],lwd=1.5)

lines(modified_boin12$ncohort*3,modified_boin12$bd.sel,type="l",pch=ltys[9],col=cs[9],lwd=1.5,lty=liness[9])
points(modified_boin12$ncohort*3,modified_boin12$bd.sel,type="p",pch=ltys[9],col=cs[9],lwd=1.5)


legend("topleft", legend=dname, col=cs, lty = liness, cex=1,lwd=rep(1.5,9),pch=ltys,bty="n",ncol=4)


dev.off()





## Figure 4
# setEPS()
# postscript(paste0(PATH,"/results/Figure4_",UTYPE,".eps"))
png(paste0(PATH,"/results/Figure4_",UTYPE,".png"), width = 7, height = 7, units = "in", res = 300)


xx = c(utpi$od.sel,
printe$od.sel,
boin$od.sel,
tepi$od.sel,
efft$od.sel,
boin12$od.sel,
modified_boin12$od.sel,
stein$od.sel,
ji3$od.sel)



par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)
plot(utpi$ncohort*3,utpi$od.sel,ylab="Favorable dose selection %",xlab="Number of patients",type="l",pch=ltys[1],col=cs[1],xlim=c(24,90),ylim=c(min(xx)*0.95,max(xx)*1.05),lwd=1.5,lty=liness[1],main=titlep)
points(utpi$ncohort*3,utpi$od.sel,type="p",pch=ltys[1],col=cs[1],lwd=1.5)

lines(printe$ncohort*3,printe$od.sel,type="l",pch=ltys[2],col=cs[2],lwd=1.5,lty=liness[2])
points(printe$ncohort*3,printe$od.sel,type="p",pch=ltys[2],col=cs[2],lwd=1.5)

lines(boin$ncohort*3,boin$od.sel,type="l",pch=ltys[3],col=cs[3],lwd=1.5,lty=liness[3])
points(boin$ncohort*3,boin$od.sel,type="p",pch=ltys[3],col=cs[3],lwd=1.5)

lines(tepi$ncohort*3,tepi$od.sel,type="l",pch=ltys[4],col=cs[4],lwd=1.5,lty=liness[4])
points(tepi$ncohort*3,tepi$od.sel,type="p",pch=ltys[4],col=cs[4],lwd=1.5)

lines(ji3$ncohort*3,ji3$od.sel,type="l",pch=ltys[5],col=cs[5],lwd=1.5,lty=liness[5])
points(ji3$ncohort*3,ji3$od.sel,type="p",pch=ltys[5],col=cs[5],lwd=1.5)

lines(boin12$ncohort*3,boin12$od.sel,type="l",pch=ltys[6],col=cs[6],lwd=1.5,lty=liness[6])
points(boin12$ncohort*3,boin12$od.sel,type="p",pch=ltys[6],col=cs[6],lwd=1.5)

lines(stein$ncohort*3,stein$od.sel,type="l",pch=ltys[7],col=cs[7],lwd=1.5,lty=liness[7])
points(stein$ncohort*3,stein$od.sel,type="p",pch=ltys[7],col=cs[7],lwd=1.5)

lines(efft$ncohort*3,efft$od.sel,type="l",pch=ltys[8],col=cs[8],lwd=1.5,lty=liness[8])
points(efft$ncohort*3,efft$od.sel,type="p",pch=ltys[8],col=cs[8],lwd=1.5)

lines(modified_boin12$ncohort*3,modified_boin12$od.sel,type="l",pch=ltys[9],col=cs[9],lwd=1.5,lty=liness[9])
points(modified_boin12$ncohort*3,modified_boin12$od.sel,type="p",pch=ltys[9],col=cs[9],lwd=1.5)



legend("topleft", legend=dname, col=cs, lty = liness, cex=1,lwd=rep(1.5,9),pch=ltys,bty="n",ncol=4)


dev.off()






## Figure 5
# setEPS()
# postscript(paste0(PATH,"/results/Figure5_",UTYPE,".eps"))
png(paste0(PATH,"/results/Figure5_",UTYPE,".png"), width = 7, height = 7, units = "in", res = 300)

xx = c(utpi$bd.pts,
printe$bd.pts,
boin$bd.pts,
tepi$bd.pts,
efft$bd.pts,
boin12$bd.pts,
modified_boin12$bd.pts,
stein$bd.pts,
ji3$bd.pts)



par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)
plot(utpi$ncohort*3,utpi$bd.pts,ylab="% patients at best dose",xlab="Number of patients",type="l",pch=ltys[1],col=cs[1],xlim=c(24,90),ylim=c(min(xx)*0.95,max(xx)*1.05),lwd=1.5,lty=liness[1],main=titlep)
points(utpi$ncohort*3,utpi$bd.pts,type="p",pch=ltys[1],col=cs[1],lwd=1.5)

lines(printe$ncohort*3,printe$bd.pts,type="l",pch=ltys[2],col=cs[2],lwd=1.5,lty=liness[2])
points(printe$ncohort*3,printe$bd.pts,type="p",pch=ltys[2],col=cs[2],lwd=1.5)

lines(boin$ncohort*3,boin$bd.pts,type="l",pch=ltys[3],col=cs[3],lwd=1.5,lty=liness[3])
points(boin$ncohort*3,boin$bd.pts,type="p",pch=ltys[3],col=cs[3],lwd=1.5)

lines(tepi$ncohort*3,tepi$bd.pts,type="l",pch=ltys[4],col=cs[4],lwd=1.5,lty=liness[4])
points(tepi$ncohort*3,tepi$bd.pts,type="p",pch=ltys[4],col=cs[4],lwd=1.5)

lines(ji3$ncohort*3,ji3$bd.pts,type="l",pch=ltys[5],col=cs[5],lwd=1.5,lty=liness[5])
points(ji3$ncohort*3,ji3$bd.pts,type="p",pch=ltys[5],col=cs[5],lwd=1.5)

lines(boin12$ncohort*3,boin12$bd.pts,type="l",pch=ltys[6],col=cs[6],lwd=1.5,lty=liness[6])
points(boin12$ncohort*3,boin12$bd.pts,type="p",pch=ltys[6],col=cs[6],lwd=1.5)

lines(stein$ncohort*3,stein$bd.pts,type="l",pch=ltys[7],col=cs[7],lwd=1.5,lty=liness[7])
points(stein$ncohort*3,stein$bd.pts,type="p",pch=ltys[7],col=cs[7],lwd=1.5)

lines(efft$ncohort*3,efft$bd.pts,type="l",pch=ltys[8],col=cs[8],lwd=1.5,lty=liness[8])
points(efft$ncohort*3,efft$bd.pts,type="p",pch=ltys[8],col=cs[8],lwd=1.5)

lines(modified_boin12$ncohort*3,modified_boin12$bd.pts,type="l",pch=ltys[9],col=cs[9],lwd=1.5,lty=liness[9])
points(modified_boin12$ncohort*3,modified_boin12$bd.pts,type="p",pch=ltys[9],col=cs[9],lwd=1.5)


legend("topleft", legend=dname, col=cs, lty = liness, cex=1,lwd=rep(1.5,9),pch=ltys,bty="n",ncol=4)


dev.off()




## Figure 6
# setEPS()
# postscript(paste0(PATH,"/results/Figure6_",UTYPE,".eps"))
png(paste0(PATH,"/results/Figure6_",UTYPE,".png"), width = 7, height = 7, units = "in", res = 300)

xx = c(utpi$poorall,
printe$poorall,
boin$poorall,
tepi$poorall,
efft$poorall,
boin12$poorall,
modified_boin12$poorall,
stein$poorall,
ji3$poorall)



par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)
plot(utpi$ncohort*3,utpi$poorall,ylab="Poor allocation %",xlab="Number of patients",type="l",pch=ltys[1],col=cs[1],xlim=c(24,90),ylim=c(min(xx)*0.95,max(xx)*1.05),lwd=1.5,lty=liness[1],main=titlep)
points(utpi$ncohort*3,utpi$poorall,type="p",pch=ltys[1],col=cs[1],lwd=1.5)

lines(printe$ncohort*3,printe$poorall,type="l",pch=ltys[2],col=cs[2],lwd=1.5,lty=liness[2])
points(printe$ncohort*3,printe$poorall,type="p",pch=ltys[2],col=cs[2],lwd=1.5)

lines(boin$ncohort*3,boin$poorall,type="l",pch=ltys[3],col=cs[3],lwd=1.5,lty=liness[3])
points(boin$ncohort*3,boin$poorall,type="p",pch=ltys[3],col=cs[3],lwd=1.5)

lines(tepi$ncohort*3,tepi$poorall,type="l",pch=ltys[4],col=cs[4],lwd=1.5,lty=liness[4])
points(tepi$ncohort*3,tepi$poorall,type="p",pch=ltys[4],col=cs[4],lwd=1.5)

lines(ji3$ncohort*3,ji3$poorall,type="l",pch=ltys[5],col=cs[5],lwd=1.5,lty=liness[5])
points(ji3$ncohort*3,ji3$poorall,type="p",pch=ltys[5],col=cs[5],lwd=1.5)

lines(boin12$ncohort*3,boin12$poorall,type="l",pch=ltys[6],col=cs[6],lwd=1.5,lty=liness[6])
points(boin12$ncohort*3,boin12$poorall,type="p",pch=ltys[6],col=cs[6],lwd=1.5)

lines(stein$ncohort*3,stein$poorall,type="l",pch=ltys[7],col=cs[7],lwd=1.5,lty=liness[7])
points(stein$ncohort*3,stein$poorall,type="p",pch=ltys[7],col=cs[7],lwd=1.5)

lines(efft$ncohort*3,efft$poorall,type="l",pch=ltys[8],col=cs[8],lwd=1.5,lty=liness[8])
points(efft$ncohort*3,efft$poorall,type="p",pch=ltys[8],col=cs[8],lwd=1.5)

lines(modified_boin12$ncohort*3,modified_boin12$poorall,type="l",pch=ltys[9],col=cs[9],lwd=1.5,lty=liness[9])
points(modified_boin12$ncohort*3,modified_boin12$poorall,type="p",pch=ltys[9],col=cs[9],lwd=1.5)


legend("topleft", legend=dname, col=cs, lty = liness, cex=1,lwd=rep(1.5,9),pch=ltys,bty="n",ncol=4)


dev.off()





## Figure 7
# setEPS()
# postscript(paste0(PATH,"/results/Figure7_",UTYPE,".eps"))
png(paste0(PATH,"/results/Figure7_",UTYPE,".png"), width = 7, height = 7, units = "in", res = 300)

xx = c(utpi$overdose,
printe$overdose,
boin$overdose,
tepi$overdose,
efft$overdose,
boin12$overdose,
modified_boin12$overdose,
stein$overdose,
ji3$overdose)



par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)
plot(utpi$ncohort*3,utpi$overdose,ylab="% patients at overdoses",xlab="Number of patients",type="l",pch=ltys[1],col=cs[1],xlim=c(24,90),ylim=c(min(xx)*0.95,max(xx)*1.05),lwd=1.5,lty=liness[1],main=titlep)
points(utpi$ncohort*3,utpi$overdose,type="p",pch=ltys[1],col=cs[1],lwd=1.5)

lines(printe$ncohort*3,printe$overdose,type="l",pch=ltys[2],col=cs[2],lwd=1.5,lty=liness[2])
points(printe$ncohort*3,printe$overdose,type="p",pch=ltys[2],col=cs[2],lwd=1.5)

lines(boin$ncohort*3,boin$overdose,type="l",pch=ltys[3],col=cs[3],lwd=1.5,lty=liness[3])
points(boin$ncohort*3,boin$overdose,type="p",pch=ltys[3],col=cs[3],lwd=1.5)

lines(tepi$ncohort*3,tepi$overdose,type="l",pch=ltys[4],col=cs[4],lwd=1.5,lty=liness[4])
points(tepi$ncohort*3,tepi$overdose,type="p",pch=ltys[4],col=cs[4],lwd=1.5)

lines(ji3$ncohort*3,ji3$overdose,type="l",pch=ltys[5],col=cs[5],lwd=1.5,lty=liness[5])
points(ji3$ncohort*3,ji3$overdose,type="p",pch=ltys[5],col=cs[5],lwd=1.5)

lines(boin12$ncohort*3,boin12$overdose,type="l",pch=ltys[6],col=cs[6],lwd=1.5,lty=liness[6])
points(boin12$ncohort*3,boin12$overdose,type="p",pch=ltys[6],col=cs[6],lwd=1.5)

lines(stein$ncohort*3,stein$overdose,type="l",pch=ltys[7],col=cs[7],lwd=1.5,lty=liness[7])
points(stein$ncohort*3,stein$overdose,type="p",pch=ltys[7],col=cs[7],lwd=1.5)

lines(efft$ncohort*3,efft$overdose,type="l",pch=ltys[8],col=cs[8],lwd=1.5,lty=liness[8])
points(efft$ncohort*3,efft$overdose,type="p",pch=ltys[8],col=cs[8],lwd=1.5)

lines(modified_boin12$ncohort*3,modified_boin12$overdose,type="l",pch=ltys[9],col=cs[9],lwd=1.5,lty=liness[9])
points(modified_boin12$ncohort*3,modified_boin12$overdose,type="p",pch=ltys[9],col=cs[9],lwd=1.5)



legend("topleft", legend=dname, col=cs, lty = liness, cex=1,lwd=rep(1.5,9),pch=ltys,bty="n",ncol=4)


dev.off()




## Figure 8
# setEPS()
# postscript(paste0(PATH,"/results/Figure8_",UTYPE,".eps"))
png(paste0(PATH,"/results/Figure8_",UTYPE,".png"), width = 7, height = 7, units = "in", res = 300)


xx = c(utpi$ov.sel,
printe$ov.sel,
boin$ov.sel,
tepi$ov.sel,
efft$ov.sel,
boin12$ov.sel,
modified_boin12$ov.sel,
stein$ov.sel,
ji3$ov.sel)

par(cex.axis=1.5, cex.lab=1.5, cex.main=1.5, cex.sub=1.5)
plot(utpi$ncohort*3,utpi$ov.sel,ylab="Overdose selection %",xlab="Number of patients",type="l",pch=ltys[1],col=cs[1],xlim=c(24,90),ylim=c(min(xx)*0.95,max(xx)*1.05),lwd=1.5,lty=liness[1],main=titlep)
points(utpi$ncohort*3,utpi$ov.sel,type="p",pch=ltys[1],col=cs[1],lwd=1.5)

lines(printe$ncohort*3,printe$ov.sel,type="l",pch=ltys[2],col=cs[2],lwd=1.5,lty=liness[2])
points(printe$ncohort*3,printe$ov.sel,type="p",pch=ltys[2],col=cs[2],lwd=1.5)

lines(boin$ncohort*3,boin$ov.sel,type="l",pch=ltys[3],col=cs[3],lwd=1.5,lty=liness[3])
points(boin$ncohort*3,boin$ov.sel,type="p",pch=ltys[3],col=cs[3],lwd=1.5)

lines(tepi$ncohort*3,tepi$ov.sel,type="l",pch=ltys[4],col=cs[4],lwd=1.5,lty=liness[4])
points(tepi$ncohort*3,tepi$ov.sel,type="p",pch=ltys[4],col=cs[4],lwd=1.5)

lines(ji3$ncohort*3,ji3$ov.sel,type="l",pch=ltys[5],col=cs[5],lwd=1.5,lty=liness[5])
points(ji3$ncohort*3,ji3$ov.sel,type="p",pch=ltys[5],col=cs[5],lwd=1.5)

lines(boin12$ncohort*3,boin12$ov.sel,type="l",pch=ltys[6],col=cs[6],lwd=1.5,lty=liness[6])
points(boin12$ncohort*3,boin12$ov.sel,type="p",pch=ltys[6],col=cs[6],lwd=1.5)

lines(stein$ncohort*3,stein$ov.sel,type="l",pch=ltys[7],col=cs[7],lwd=1.5,lty=liness[7])
points(stein$ncohort*3,stein$ov.sel,type="p",pch=ltys[7],col=cs[7],lwd=1.5)

lines(efft$ncohort*3,efft$ov.sel,type="l",pch=ltys[8],col=cs[8],lwd=1.5,lty=liness[8])
points(efft$ncohort*3,efft$ov.sel,type="p",pch=ltys[8],col=cs[8],lwd=1.5)

lines(modified_boin12$ncohort*3,modified_boin12$ov.sel,type="l",pch=ltys[9],col=cs[9],lwd=1.5,lty=liness[9])
points(modified_boin12$ncohort*3,modified_boin12$ov.sel,type="p",pch=ltys[9],col=cs[9],lwd=1.5)



legend("topleft", legend=dname, col=cs, lty = liness, cex=1,lwd=rep(1.5,9),pch=ltys,bty="n",ncol=4)


dev.off()



}



