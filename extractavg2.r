library(data.table)



extractavg<-function(path){
for(i in 1:NDOSE){
utpi= fread(paste0(i,path))
KK = ncol(utpi)-1
if (i==1){
utpix = utpi[,1:KK]
}else{
utpix = utpi[,1:KK] + utpix
}
}
utpix = utpix / NDOSE
return(cbind(utpix,utpi$design))
}

utpi= extractavg("\\utpi_random.csv")
boin=extractavg("\\BOINET_random.csv")
efft=temp("\\efftox_random.csv")
tepi=extractavg("\\tepi_random.csv")
boin12 = extractavg("\\boin12_random.csv")
ji3 = extractavg("\\ji3_random.csv")
printe = extractavg("\\printe_random.csv")
stein = extractavg("\\stein_random.csv")

tblall = NULL
tblall = rbind(tblall, utpi)
tblall = rbind(tblall, boin)
tblall = rbind(tblall, efft)
tblall = rbind(tblall, tepi)
tblall = rbind(tblall, boin12)
tblall = rbind(tblall, ji3)
tblall = rbind(tblall, printe)
tblall = rbind(tblall, stein)


tblall=as.data.frame(tblall)
tblall = tblall[which(tblall$ncohort %in% SSIZERANGE),]
fwrite(tblall,"simres.csv")

