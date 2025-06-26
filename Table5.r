
library(data.table)
library(tidyverse)

ub.setting <- UB

tbl = fread(paste0(PATH,'/intermediate/NDOSE10/simres', ub.setting, '.csv'))

tbl$Design = tbl$V2
tbl$Design <- ifelse(tbl$Design == "utpi", "uTPI", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "boinet", "BOIN-ET", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "ji3", "Joint3+3", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "stein", "STEIN", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "printe", "PRINTE", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "tepi", "TEPI", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "boin12", "BOIN12", tbl$Design)
tbl$Design <- ifelse(tbl$Design == paste0("modified", ub.setting, "_boin12"), "M-BOIN12", tbl$Design)
tbl$Design <- ifelse(tbl$Design == "efftox", "EffTox", tbl$Design)
		   

p1 = tbl %>% filter(rtype ==1 & utype ==1)
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
p1 = p1 %>% select('Sample Size',Design, "OBD Sel%","FD Sel%","OBD Pts%","Poor Pts%","OV Pts%","OV Sel%",CompRank)




p2 = tbl %>% filter(rtype ==1 & utype ==2)

p2=as.data.table(p2)
preferred.order = c("uTPI","BOIN-ET","Joint3+3","STEIN","PRINTE", "TEPI","BOIN12", "M-BOIN12")
p2 = p2[preferred.order, on="Design"]
p2 = p2 %>% arrange(ncohort)


p2$temp = paste0(p2$ncohort,p2$Design)
p2[,bdselr := rank(bd.sel * (-1)),by=ncohort]
p2[,odselr := rank(od.sel * (-1)),by=ncohort]
p2[,bdptsr := rank(bd.pts * (-1)),by=ncohort]
p2[,poorallr := rank(poorall ),by=ncohort]
p2[,overdoser := rank(overdose),by=ncohort]
p2[,ovselr := rank(ov.sel),by=ncohort]

p2[,avgrank := (bdselr+odselr+bdptsr+poorallr+overdoser+ovselr)/6]
p2[,CompRank := rank(bdselr+odselr+bdptsr+poorallr+overdoser+ovselr, ties.method ="min"),by=ncohort]
p2$avgrank = format(round(p2$avgrank, 2), nsmall = 2)
p2$CompRank = paste0(p2$CompRank,"(",p2$avgrank,")")

p2$bd.sel = format(round(p2$bd.sel, 2), nsmall = 2)
p2$od.sel = format(round(p2$od.sel, 2), nsmall = 2)
p2$bd.pts = format(round(p2$bd.pts, 2), nsmall = 2)
p2$poorall = format(round(p2$poorall, 2), nsmall = 2)
p2$overdose = format(round(p2$overdose, 2), nsmall = 2)
p2$ov.sel = format(round(p2$ov.sel, 2), nsmall = 2)

p2$'Sample Size' = p2$ncohort * 3
p2 = p2 %>% rename("OBD Sel%" = "bd.sel","FD Sel%"="od.sel",	"OBD Pts%"="bd.pts",	"Poor Pts%"="poorall",	"OV Pts%"="overdose",	"OV Sel%"="ov.sel")
p2 = p2 %>% select('Sample Size',Design, "OBD Sel%","FD Sel%","OBD Pts%","Poor Pts%","OV Pts%","OV Sel%",CompRank)


fwrite((rbind(p2,p1)),paste0("results/", ub.setting, "/Table5.csv"))