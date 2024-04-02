
Source code for manuscript "Comparative Review of Novel Model-Assisted Designs for Phase I/II Clinical Trials"

Haolun Shi, Ruitao Lin, Xiaolei Lin


# Introduction

The code has been written using R-3.6.3 (platform: x86_64-w64-mingw32
x64 (64-bit) with package versions data.table_1.14.2 and Iso_0.0-18.1.

To reproduce the tables in the manuscript, run the main_part1.r then the main_part2.r. The script
will produce all the intermediate results (in csv format) and Figures. A
detailed description of each R code are given as follows.

-   main_part1.r: The script will produce all the intermediate results (in csv
    format) and Figures.

-   main_part2.r: The script will produce all the numerical tables (in csv
    format) and Figures based on the intermediate results. 

-   Figure1.r: The script will reproduce Figure 1 in the manuscript.

-   Figure2.r: The script will reproduce Figure 2 in the manuscript.

-   Figure3-8.r: The script will reproduce Figures 3-8 in the
    manuscript.

-   FigureS1-S10.r: The script will reproduce Figures in the
    supplementary material.

-   Table2.r: The script will reproduce Table 2 in the
    manuscript.
	
-   Table3.r: The script will reproduce Table 3 in the
    manuscript.

-   Table4.r: The script will reproduce Table 4 in the
    manuscript.
	
-   Table5.r: The script will reproduce Table 5 in the
    manuscript.
	
-   Table6.r: The script will reproduce Table 6 in the
    manuscript.	

-   boin12.R: The script for generating operating characteristics of the
    BOIN12 design.

-   BOINET.R: The script for generating operating characteristics of the
    BOIN-ET design.

-   efftox.r: The script for generating operating characteristics of the
    EffTox design.

-   Ji3.R: The script for generating operating characteristics of the
    Joint i3+3 design.

-   PITE.R: The script for generating operating characteristics of the
    PRINTE design.

-   STEIN.R: The script for generating operating characteristics of the
    STEIN design.

-   TEPI.R: The script for generating operating characteristics of the
    TEPI design.

-   uTPI.R: The script for generating operating characteristics of the
    uTPI design.

-   extractavg.r: The script for extracting and averaging the operating
    characteristics of the model-assisted design into a single csv files
    named "simres.csv".

-   extractavg2.r: The script for extracting and averaging the operating
    characteristics of the model-assisted design and EffTox design into
    a single csv files named "simres.csv".

-   extractavg_sens.r: The script for extracting and averaging the
    operating characteristics of the model-assisted design into a single
    csv files named "simres_sensitivity.csv" for sensitivity analyses.

Note that the running time is long and we have stored all the
intermediate results under the folder "intermediate" .
Tables 1 and 2 are text-based tables and contain no numbers, thus are
omitted.


Description of the folder structures are as follows. 

-   results/: The folder contains all the numerical tables and figures 

-   intermediate/: The folder contains all the intermediate results for simulation    

-   intermediate/NDOSE3/: The folder contains the intermediate results for simulation
    with the number of dose level $J=3$.

-   intermediate/NDOSE5/: The folder contains the intermediate results for simulation
    with the number of dose level $J=5$.

-   intermediate/NDOSE5_START2/: The folder contains the intermediate results for
    simulation with the number of dose level $J=5$ and starting dose
    level 2.

-   intermediate/NDOSE10/: The folder contains the intermediate results for
    simulation with the number of dose level $J=10$.

-   intermediate/sensitivity/: The folder contains the intermediate results for
    simulation on sensitivity analyses. The subfolder within corresponds
    to the simulation with a certain target efficacy rate $\zeta_E$ and
    target toxicity rate $\phi_T$ (subfolder name is simply the
    $\zeta_E$ and $\phi_T$ connected by a "\_".

Each folder contains a number of subfolders indexed from 1 to $J$, which
correspond to the scenario where the OBD is fixed at the index. For
example, the subfolder NDOSE3/2 corresponds to the simulation scenario
where the OBD is fixed at the second dose level in the simulation where
there are $J=3$ dose levels.

For questions, comments or remarks about the code, please contact the
first author of the paper.

# Additional Guidance for User-Specified Scenarios

In practice, a user may prefer to specify parameters that are different
from the ones in the manuscript. Moreover, the user may prefer to set up
fixed scenarios for comparison. In this section, we provide additional
guidance on how to compare the user-specified designs and setting up
user-specified ground truth.

Note that this section only serve as an extra guidance note for readers
that are interested to use the code to carry on their own simulation
studies. This section does not correspond to or reproduce any of the
results in the manuscript. The main.r is already suffices to reproduce
all tables and figures.

## Setting Up Ground Truth: Random Scenarios

In the case of user-specified random scenarios, the user may simply set
the global variables `TARGET.E` and `TARGET.T` before running the
simulations. For example, if a user prefer to generate random scenarios
with $\zeta_T = 0.25$ and $\phi_T = 0.35$ and evaluate the operating
characteristics of BOIN12 design with 5 dose levels and number of
cohorts being 30, the following code can be used.

``` {.r language="R"}
PATH = getwd()

NDOSE <<- 5
dir.create(paste0(PATH,'/NDOSE5/'), showWarnings = FALSE)
setwd(paste0(PATH,'/NDOSE5/'))

for(i in 1:NDOSE){
dir.create(paste0(i,'/'), showWarnings = FALSE)
}

STARTD <<- 1
TARGET.E <<- 0.25
TARGET.T <<- 0.35
SSIZERANGE <<- c(30)
source(paste0(PATH,'/boin12.R'))
```

## Setting Up Ground Truth: Fixed Scenarios

In the case of user-specified fixed scenarios, a fixed probability
vectors needs to be supplied for the toxicity probability and efficacy
probability. To do this, the user can simply replace the `simprob()`
function within the source codes of the designs.

For example, if a user prefer to evaluate operating characteristics
where the toxicity and efficacy probabilities for the five doses are
$(0.1,0.2,0.3,0.4,0.4)$ and $(0.4,0.5,0.6,0.6,0.6)$, the following code
can be used to replace the `simprob()` function. of cohorts being 30,
the following code can be used to replace the `simprob()` function
within the source codes of the designs.

``` {.r language="R"}
simprob<-function(ndose,targetE,targetT,u1,u2 ,randomtype){
  pT = c(0.1,0.2,0.3,0.4,0.4)
  pE = c(0.4,0.5,0.6,0.6,0.6)
  obd = 3
  mtd = 2
  return(list(pE = pE,pT=pT,obd=obd,mtd=mtd))
}
```

## Setting Design Parameters

We describe how to set up the design parameters in the following
subsections. For BOIN-ET and EffTox designs, many of the key parameters
are calculated from external software/codes. The code for BOIN-ET can be
obtained upon request from Takeda, K., Taguri, M., and Morita, S, who
are the authors of the paper. As for the software for implementing the
EffTox design, we refer the readers to the website
`https://biostatistics.mdanderson.org/SoftwareDownload/`.

### BOIN12 Design Parameters

We describe the design parameters that underlies the arguments for the
function `get.oc.obd()` in the BOIN12.R as follows.

-   `targetT`: target toxicity probability, i.e., the $\phi_T$ in the
    manuscript.

-   `targetE`: lower limit for the efficacy rate, i.e., the $\zeta_E$ in
    the manuscript.

-   `ncohort`: the number of cohorts.

-   `cohortsize`: the size of a cohort.

-   `startdose`: the starting dose level.

-   `cutoff.eli.T`: the posterior probability cutoff for toxicity dose
    elimination rule.

-   `cutoff.eli.E`: the posterior probability cutoff for futility dose
    elimination rule.

-   `u1`: the utility parameter $w_{11}$, in the scale of 0 to 100.

-   `u2`: the utility parameter $w_{00}$, in the scale of 0 to 100.

-   `ntrial`: the number of random trial replications.

-   `utilitytype`: a overriding argument for controlling the type of
    utility parameters. If set as 1, then $(w_{11},w_{00}) = (0.6,0.4)$.
    If set as 2, $(w_{11},w_{00}) = (1,0)$. Otherwise, the
    user-specified values for `u1` and `u2` are used.

The early stage sample size cutoff is set as the recommended default
values, and rarely needs to be modified. The rest of the design
parameters are intermediate and depend upon the user-specified
parameters listed above.

If the user wish to specify utility structure, he may simply modify the
values in the following section within the source code and set
`utilitytype = 1`.

``` {.r language="R"}
if (utilitytype==1){
u1 = 60
u2 = 40
}
```

### uTPI Design Parameters

We describe the design parameters that underlies the arguments for the
function `get.oc()` in the uTPI.R as follows.

-   `targetT`: target toxicity probability, i.e., the $\phi_T$ in the
    manuscript.

-   `targetE`: lower limit for the efficacy rate, i.e., the $\zeta_E$ in
    the manuscript.

-   `ncohort`: the number of cohorts.

-   `cohortsize`: the size of a cohort.

-   `startdose`: the starting dose level.

-   `cutoff.eli.T`: the posterior probability cutoff for toxicity dose
    elimination rule.

-   `cutoff.eli.E`: the posterior probability cutoff for futility dose
    elimination rule.

-   `ntrial`: the number of random trial replications.

-   `utilitytype`: a overriding argument for controlling the type of
    utility parameters. If set as 1, then $(w_{11},w_{00}) = (0.6,0.4)$.
    If set as 2, $(w_{11},w_{00}) = (1,0)$.

The early stage sample size cutoff is set as the recommended default
values, and rarely needs to be modified. The rest of the design
parameters are intermediate and depend upon the user-specified
parameters listed above.

If the user wish to specify utility structure, he may simply modify the
values in the following section within the source code and set
`utilitytype = 1`.

``` {.r language="R"}
if (utilitytype==1){
u1 = 60
u2 = 40
}
```

### TEPI Design Parameters

We describe the design parameters that underlies the arguments for the
function `get.oc()` in the TEPI.R as follows.

-   `targetT`: target toxicity probability, i.e., the $\phi_T$ in the
    manuscript.

-   `targetE`: lower limit for the efficacy rate, i.e., the $\zeta_E$ in
    the manuscript.

-   `ncohort`: the number of cohorts.

-   `cohortsize`: the size of a cohort.

-   `startdose`: the starting dose level.

-   `cutoff.eli.T`: the posterior probability cutoff for toxicity dose
    elimination rule.

-   `cutoff.eli.E`: the posterior probability cutoff for futility dose
    elimination rule.

-   `ntrial`: the number of random trial replications.

-   `utilitytype`: a overriding argument for controlling the type of
    utility parameters. If set as 1, then $(w_{11},w_{00}) = (0.6,0.4)$.
    If set as 2, $(w_{11},w_{00}) = (1,0)$.

Note that the TEPI design requires a mapping table for dose assignment
decision. These can be modified and specified in line 220 - 223 as
follows (here, we set the table to be hinge upon the target efficacy and
toxicity probability).

``` {.r language="R"}
effint_l<-c(0,TARGET.E,TARGET.E+0.2,TARGET.E+0.4)
effint_u<-c(TARGET.E,TARGET.E+0.2,TARGET.E+0.4,1)
toxint_l<-c(0,0.15,TARGET.T,TARGET.T+0.05)
toxint_u<-c(0.15,TARGET.T,TARGET.T+0.05,1)
```

The 4 variables in the code above set up the boundaries within the
mapping tables. Setting `TARGET.E` to be 0.2 and `TARGET.T` to be 0.3
reproduce Table 1 in the manuscript.

If the user wish to specify utility structure, he may simply modify the
values in the following section within the source code and set
`utilitytype = 1`.

``` {.r language="R"}
if (utilitytype==1){
u1 = 60
u2 = 40
}
```

### PRINTE Design Parameters

We describe the design parameters that underlies the arguments for the
function `PITE_sim()` in the PITE.R as follows.

-   `pT`: target toxicity probability, i.e., the $\phi_T$ in the
    manuscript.

-   `pE`: the target efficacy probability, i.e., the $\phi_E$ in the
    manuscript.

-   `zetaE`: lower limit for the efficacy rate, i.e., the $\zeta_E$ in
    the manuscript.

-   `ncohort`: the number of cohorts.

-   `cohortsize`: the size of a cohort.

-   `startdose`: the starting dose level.

-   `eps1`: Width of the subrectangle $\epsilon$.

-   `eps2`: Width of the subrectangle $\epsilon$.

-   `psafe`: the posterior probability cutoff for toxicity dose
    elimination rule.

-   `pfutility`: the posterior probability cutoff for futility dose
    elimination rule.

-   `nsimul`: the number of random trial replications.

-   `utilitytype`: a overriding argument for controlling the type of
    utility parameters. If set as 1, then $(w_{11},w_{00}) = (0.6,0.4)$.
    If set as 2, $(w_{11},w_{00}) = (1,0)$.

If the user wish to specify utility structure, he may simply modify the
values in the following section within the source code and set
`utilitytype = 1`.

``` {.r language="R"}
if (utilitytype==1){
u1 = 60
u2 = 40
}
```

### Joint i3+3 Design Parameters

We describe the design parameters that underlies the arguments for the
function `PITE_sim()` in the PITE.R as follows.

-   `pT`: target toxicity probability, i.e., the $\phi_T$ in the
    manuscript.

-   `pE`: the target efficacy probability, i.e., the $\phi_E$ in the
    manuscript.

-   `zetaE`: lower limit for the efficacy rate, i.e., the $\zeta_E$ in
    the manuscript.

-   `ncohort`: the number of cohorts.

-   `cohortsize`: the size of a cohort.

-   `startdose`: the starting dose level.

-   `eps1`: Width of the subrectangle $\epsilon$.

-   `eps2`: Width of the subrectangle $\epsilon$.

-   `psafe`: the posterior probability cutoff for toxicity dose
    elimination rule.

-   `pfutility`: the posterior probability cutoff for futility dose
    elimination rule.

-   `nsimul`: the number of random trial replications.

-   `utilitytype`: a overriding argument for controlling the type of
    utility parameters. If set as 1, then $(w_{11},w_{00}) = (0.6,0.4)$.
    If set as 2, $(w_{11},w_{00}) = (1,0)$.

If the user wish to specify utility structure, he may simply modify the
values in the following section within the source code and set
`utilitytype = 1`.

``` {.r language="R"}
if (utilitytype==1){
u1 = 60
u2 = 40
}
```

### STEIN Design Parameters

We describe the design parameters that underlies the arguments for the
function `get.oc()` in the STEIN.R as follows.

-   `targetT`: target toxicity probability, i.e., the $\phi_T$ in the
    manuscript.

-   `targetE`: lower limit for the efficacy rate, i.e., the $\zeta_E$ in
    the manuscript.

-   `psi1`: the highest efficacy probability deemed as inefficacious,
    i.e., the $\psi_{1,E}$ in the manuscript.

-   `psi2`: the lowest efficacy probability deemed as highly promising,
    i.e., the $\psi_{2,E}$ in the manuscript.

-   `ncohort`: the number of cohorts.

-   `cohortsize`: the size of a cohort.

-   `startdose`: the starting dose level.

-   `cutoff.eli.T`: the posterior probability cutoff for toxicity dose
    elimination rule.

-   `cutoff.eli.E`: the posterior probability cutoff for futility dose
    elimination rule.

-   `ntrial`: the number of random trial replications.

-   `utilitytype`: a overriding argument for controlling the type of
    utility parameters. If set as 1, then $(w_{11},w_{00}) = (0.6,0.4)$.
    If set as 2, $(w_{11},w_{00}) = (1,0)$.

If the user wish to specify utility structure, he may simply modify the
values in the following section within the source code and set
`utilitytype = 1`.

``` {.r language="R"}
if (utilitytype==1){
u1 = 60
u2 = 40
}
```

The value of the design parameter $\eta_1$ is calculated from `psi1` and
`psi2` at line 214 in the source code as follows. Note that in a
comparison studies, this value needs to be aligned to be close to the
$\phi_E$ in the joint i3+3 and PRINTE design (or `pE` as the functional
argument).

``` {.r language="R"}
psi<-log((1-psi1)/(1-psi2))/log(psi2*(1-psi1)/psi1/(1-psi2))
```
