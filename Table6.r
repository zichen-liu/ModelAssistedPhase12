
library(data.table)
library(tidyverse)



tbl = fread(paste0(PATH,'/intermediate/sensitivity/simres_sensitivity.csv'))

tbl$Design = tbl$V2
tbl$Design <- ifelse(tbl$Design == "utpi", "uTPI", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "boinet", "BOIN-ET", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "ji3", "Joint3+3", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "stein", "STEIN", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "printe", "PRINTE", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "tepi", "TEPI", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "boin12", "BOIN12", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "modified_boin12", "M-BOIN12", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "efftox", "EffTox", tbl$Design)
		   

p1 = tbl %>% filter(rtype ==1 & utype ==1 & ncohort == 18 & toxtarget == 0.3 & efftarget == 0.1)
p2 = tbl %>% filter(rtype ==1 & utype ==1 & ncohort == 18 & toxtarget == 0.3 & efftarget == 0.2)
p3 = tbl %>% filter(rtype ==1 & utype ==1 & ncohort == 18 & toxtarget == 0.3 & efftarget == 0.3)
p4 = tbl %>% filter(rtype ==1 & utype ==1 & ncohort == 18 & toxtarget == 0.3 & efftarget == 0.4)
p5 = tbl %>% filter(rtype ==1 & utype ==1 & ncohort == 18 & toxtarget == 0.2 & efftarget == 0.2)
p6 = tbl %>% filter(rtype ==1 & utype ==1 & ncohort == 18 & toxtarget == 0.4 & efftarget == 0.2)

tableformat<-function(p1,toxtarget,efftarget){

p1=as.data.table(p1)
preferred.order = c("uTPI","BOIN-ET","Joint3+3","STEIN","PRINTE", "TEPI","BOIN12", "M-BOIN12")
p1 = p1[preferred.order, on="Design"]
p1 = p1 %>% arrange(ncohort)


p1$temp = paste0(p1$ncohort,p1$Design)
p1[,bdselr := rank(bd.sel * (-1)),by=ncohort]
p1[,odselr := rank(od.sel * (-1)),by=ncohort]
p1[,bdptsr := rank(bd.pts * (-1)),by=ncohort]
p1[,poorallr := rank(poorall ),by=ncohort]
p1[,overdoser := rank(overdose),by=ncohort]
p1[,ovselr := rank(ov.sel),by=ncohort]

p1[,avgrank := (bdselr+odselr+bdptsr+poorallr+overdoser+ovselr)/6]
p1[,CompRank := rank(bdselr+odselr+bdptsr+poorallr+overdoser+ovselr, ties.method ="min"),by=ncohort]
p1$avgrank = format(round(p1$avgrank, 2), nsmall = 2)
p1$CompRank = paste0(p1$CompRank,"(",p1$avgrank,")")

p1$bd.sel = format(round(p1$bd.sel, 2), nsmall = 2)
p1$od.sel = format(round(p1$od.sel, 2), nsmall = 2)
p1$bd.pts = format(round(p1$bd.pts, 2), nsmall = 2)
p1$poorall = format(round(p1$poorall, 2), nsmall = 2)
p1$overdose = format(round(p1$overdose, 2), nsmall = 2)
p1$ov.sel = format(round(p1$ov.sel, 2), nsmall = 2)

p1$'Sample Size' = p1$ncohort * 3
p1 = p1 %>% rename("OBD Sel%" = "bd.sel","FD Sel%"="od.sel",	"OBD Pts%"="bd.pts",	"Poor Pts%"="poorall",	"OV Pts%"="overdose",	"OV Sel%"="ov.sel")
p1$'phi_T' = p1$toxtarget
p1$'zeta_E' = p1$efftarget

p1 = p1 %>% select('phi_T','zeta_E',Design, "OBD Sel%","FD Sel%","OBD Pts%","Poor Pts%","OV Pts%","OV Sel%",CompRank)

return(p1)
}


output = rbind(
tableformat(p1,0.3,0.1),
tableformat(p2,0.3,0.2),
tableformat(p3,0.3,0.3),
tableformat(p4,0.3,0.4),
tableformat(p5,0.2,0.2),
tableformat(p6,0.4,0.2)
)

fwrite(output,paste0("results/Table6.csv"))

