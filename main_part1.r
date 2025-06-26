# The following code reproduce the intermediate results of the paper "Comparative Review of Novel Model-Assisted Designs for Phase I/II Clinical Trials" by Haolun Shi, Ruitao Lin, and Xiaolei Lin. 
# The code needs to be run first in order to produce all the figures and tables in the manuscript

PATH = getwd()

# Benchmark utility setting for modified BOIN12
# Select which modified benchmark utility to use instead of
# the standard BOIN12 ub = uu + (100 - uu) / 2
# 1: ub = uu
# 2: ub = uu + (100 - uu) / 4
# 3: ub = 100 - (100 - uu) / 4
UB <- 2

# Run Number of dose = 5, Starting dose = 1
NDOSE <<- 5
# dir.create(paste0(PATH,'/intermediate/NDOSE5/'), showWarnings = FALSE)
setwd(paste0(PATH,'/intermediate/NDOSE5/'))

# for(i in 1:NDOSE){
# dir.create(paste0(i,'/'), showWarnings = FALSE)
# }

STARTD <<- 1
TARGET.E <<- 0.2 
TARGET.T <<- 0.3 
SSIZERANGE <<- seq(8,30,by=2)
# source(paste0(PATH,'/boin12.R'))
source(paste0(PATH,'/modified_boin12.R'))
# source(paste0(PATH,'/BOINET.R'))
# source(paste0(PATH,'/Ji3.R'))
# source(paste0(PATH,'/PITE.R'))
# source(paste0(PATH,'/STEIN.R'))
# source(paste0(PATH,'/TEPI.R'))
# source(paste0(PATH,'/uTPI.R'))
# source(paste0(PATH,'/efftox.r'))
# setwd(paste0(PATH,'/NDOSE5/'))
source(paste0(PATH,'/extractavg2.r'))



# Run Number of dose = 5, Starting dose = 2
NDOSE <<- 5
# dir.create(paste0(PATH,'/intermediate/NDOSE5_START2/'), showWarnings = FALSE)
setwd(paste0(PATH,'/intermediate/NDOSE5_START2/'))

# for(i in 1:NDOSE){
# dir.create(paste0(i,'/'), showWarnings = FALSE)
# }

STARTD <<- 2
TARGET.E <<- 0.2 
TARGET.T <<- 0.3 
SSIZERANGE <<- seq(8,30,by=2)
# source(paste0(PATH,'/boin12.R'))
source(paste0(PATH,'/modified_boin12.R'))
# source(paste0(PATH,'/BOINET.R'))
# source(paste0(PATH,'/Ji3.R'))
# source(paste0(PATH,'/PITE.R'))
# source(paste0(PATH,'/STEIN.R'))
# source(paste0(PATH,'/TEPI.R'))
# source(paste0(PATH,'/uTPI.R'))
# source(paste0(PATH,'/efftox.r'))
# setwd(paste0(PATH,'/NDOSE5_START2/'))
source(paste0(PATH,'/extractavg2.r'))


# Run Number of dose = 3, Starting dose = 1
NDOSE <<- 3
# dir.create(paste0(PATH,'/intermediate/NDOSE3/'), showWarnings = FALSE)
setwd(paste0(PATH,'/intermediate/NDOSE3/'))

# for(i in 1:NDOSE){
# dir.create(paste0(i,'/'), showWarnings = FALSE)
# }

STARTD <<- 1
TARGET.E <<- 0.2 
TARGET.T <<- 0.3 
SSIZERANGE <<- c(6,12)
# source(paste0(PATH,'/boin12.R'))
source(paste0(PATH,'/modified_boin12.R'))
# source(paste0(PATH,'/BOINET.R'))
# source(paste0(PATH,'/Ji3.R'))
# source(paste0(PATH,'/PITE.R'))
# source(paste0(PATH,'/STEIN.R'))
# source(paste0(PATH,'/TEPI.R'))
# source(paste0(PATH,'/uTPI.R'))
# setwd(paste0(PATH,'/NDOSE3/'))
source(paste0(PATH,'/extractavg.r'))






# Run Number of dose = 10, Starting dose = 1
NDOSE <<- 10
# dir.create(paste0(PATH,'/intermediate/NDOSE10/'), showWarnings = FALSE)
setwd(paste0(PATH,'/intermediate/NDOSE10/'))

# for(i in 1:NDOSE){
# dir.create(paste0(i,'/'), showWarnings = FALSE)
# }

STARTD <<- 1
TARGET.E <<- 0.2 
TARGET.T <<- 0.3 
SSIZERANGE <<- c(18,36)

# source(paste0(PATH,'/boin12.R'))
source(paste0(PATH,'/modified_boin12.R'))
# source(paste0(PATH,'/BOINET.R'))
# source(paste0(PATH,'/Ji3.R'))
# source(paste0(PATH,'/PITE.R'))
# source(paste0(PATH,'/STEIN.R'))
# source(paste0(PATH,'/TEPI.R'))
# source(paste0(PATH,'/uTPI.R'))
# setwd(paste0(PATH,'/NDOSE10/'))
source(paste0(PATH,'/extractavg.r'))





# Run - Sensitivity Analysis
NDOSE <<- 5
STARTD <<- 1
SSIZERANGE <<- c(12,18)

ERANGE = c(0.1,0.2,0.3,0.4)

TRANGE = c(0.2,0.3,0.4)

# dir.create(paste0(PATH,'/intermediate/sensitivity/'), showWarnings = FALSE)
setwd(paste0(PATH,'/intermediate/sensitivity/'))
for (ee in ERANGE){
for (tt in TRANGE){
	
	setwd(paste0(PATH,'/intermediate/sensitivity/'))
	dir.create(paste0(ee,'_',tt,'/'), showWarnings = FALSE)
	setwd(paste0(ee,'_',tt,'/'))
	for(i in 1:NDOSE){
	dir.create(paste0(i,'/'), showWarnings = FALSE)
	}

	TARGET.E <<- ee
	TARGET.T <<- tt 
	# source(paste0(PATH,'/boin12.R'))
	source(paste0(PATH,'/modified_boin12.R'))
	# source(paste0(PATH,'/BOINET.R'))
	# source(paste0(PATH,'/Ji3.R'))
	# source(paste0(PATH,'/PITE.R'))
	# source(paste0(PATH,'/STEIN.R'))
	# source(paste0(PATH,'/TEPI.R'))
	# source(paste0(PATH,'/uTPI.R'))
}
}

source(paste0(PATH,'/extractavg_sens.r'))



