# The following code reproduce the all the figures and tables of the paper "Comparative Review of Novel Model-Assisted Designs for Phase I/II Clinical Trials" by Haolun Shi, Ruitao Lin, and Xiaolei Lin. 
# The code 'main_part1.r' needs to be run first in order to run this script


PATH = getwd()
setwd(paste0(PATH,'/intermediate/NDOSE5/'))
setwd("../../")
# source(paste0(PATH,'/Figure1.r'))
# source(paste0(PATH,'/Figure2.r'))
source(paste0(PATH,'/Figure3-8.r'))

source(paste0(PATH,'/FigureS1-S10.r'))
setwd("../../..")

source(paste0(PATH,'/Table2.r'))
source(paste0(PATH,'/Table3.r'))
source(paste0(PATH,'/Table4.r'))
source(paste0(PATH,'/Table5.r'))
source(paste0(PATH,'/Table6.r'))

