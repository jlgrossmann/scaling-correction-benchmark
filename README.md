## scaling-correction-benchmark
This repository contains code and data to accompany [JL Großmann et al: "Critical evaluation of assessor difference correction approaches in sensory analysis"](https://doi.org/10.1016/j.foodqual.2022.104792). 
The repository is intented to illustrate the scaling correction, model fitting, permutation and simulation strategies employed in the aforementioned paper.

### `permutation_simulation.R`
This script runs the permutation or simulation and saves the results to disk for later analysis. 
The options/parameters for the simulation and permutation approaches are documented in the code. 
Parameters can be either hard-coded into the script (some defaults are given) or the script can be run from the command line using `Rscript` and parameters can be specified as arguments. 
Example: `Rscript permutation_simulation.R --mode=permutation --ptyp=ass --nsim=10000 --ratt 13` runs restricted permutation on the 13th sensory attribute in the dataset (see below) 10000 times. 

### `helper_functions.R`
This file contains functions for scaling correction (standardisation, ten Berge method, assessor model) and ANOVA fitting.

### `MAManalysis.R`
This script fits the balanced mixed assessor model (MAM). The script was modified from the [`SensMixed` package](https://cran.r-project.org/web/packages/SensMixed/index.html), which was published by A Kuznetsova, PB Brockhoff and RHB Christensen under a GPL-2 license. 

### `evaluation_plot.R`
This script shows an example for how to further process the results produced by `permutation_simulation.R` and obtain p-values from the permutation approach. 

### `cheese.csv`
The dataset used by [Romano et al (2008)](https://doi.org/10.1016/j.foodqual.2007.06.008) and as an example in this study. 

