
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


NDOSE <<- 5
STARTD <<- 1
SSIZERANGE <<- c(12,18)

ERANGE = c(0.1,0.2,0.3,0.4)

TRANGE = c(0.2,0.3,0.4)


setwd(paste0(PATH,'intermediate/sensitivity/'))
for (ee in ERANGE){
for (tt in TRANGE){
	
	setwd(paste0(PATH,'intermediate/sensitivity/'))
	
	setwd(paste0(ee,'_',tt,'/'))
	
	TARGET.E <<- ee
	TARGET.T <<- tt 
		
		

	utpi= extractavg("\\utpi_random.csv")
	boin=extractavg("\\BOINET_random.csv")
	tepi=extractavg("\\tepi_random.csv")
	boin12 = extractavg("\\boin12_random.csv")
	ji3 = extractavg("\\ji3_random.csv")
	printe = extractavg("\\printe_random.csv")
	stein = extractavg("\\stein_random.csv")

	tblall = NULL
	tblall = rbind(tblall, utpi)
	tblall = rbind(tblall, boin)
	tblall = rbind(tblall, tepi)
	tblall = rbind(tblall, boin12)
	tblall = rbind(tblall, ji3)
	tblall = rbind(tblall, printe)
	tblall = rbind(tblall, stein)


	tblall=as.data.frame(tblall)
	tblall = tblall[which(tblall$ncohort %in% c(12,18)),]
	fwrite(tblall,"simres.csv")

}
}








tblsens = NULL

for (ee in ERANGE){
for (tt in TRANGE){
	
	setwd(paste0(PATH,'intermediate/sensitivity/'))
	
	setwd(paste0(ee,'_',tt,'/'))
	
	tblall = fread("simres.csv")

	tblall$efftarget = ee 
	tblall$toxtarget = tt
	tblsens = rbind(tblsens,tblall)
}
}

setwd(paste0(PATH,'intermediate/sensitivity/'))
fwrite(tblsens,'simres_sensitivity.csv')
