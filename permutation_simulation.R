# supplementary R code to accompany
# JL Großmann et al: "Critical evaluation of assessor difference correction approaches in sensory analysis"
# the following code is intented to illustrate the scaling correction, model fitting, permutation and 
# simulation strategies employed in the aforementioned paper 

# helper_functions.R contains functions for scaling correction and model fitting
source('helper_functions.R')
# MAManalysis.R is derived from the SensMixed R package by A Kuznetsova, PB Brockhoff and RHB Christensen 
source('MAManalysis.R')
starttime = format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
library(getopt)  # argument parsing

# this script can be used to replicate the results of the permutation and simulation
# approaches laid out in JL Großmann et al.
# working directory must contain this script, helper_functions.R, MAManalysis.R and cheese.csv
# options can either be hard-coded below or set as arguments in command line:
# Rscript permutation_simulation.R --mode=permutation --nsim=10000 --ratt 13 --ptyp=ass
#   to perform restricted permutation on Fattiness atribute, 10000 iterations
# Rscript permutation_simulation.R --mode=simulation --nsim=10000 --simm=random_equal
#   to perform simulation with equal sigma between assessors, 10000 iterations



######## PARAMETERS #########
# strategy to employ to obtain null distributions
mode. = 'permutation'  # options: 'permutation' - will permute selected attribute of cheese dataset
                       #          'simulation' - will generate random data 
n_sim = 1e4  # Number of permutation / simulation iterations

# for permutation
perm_type = 'ass'  # how to permute 
                   # options: 'ass' - "restricted permutation"
                   #          'assrepl'
                   #          'all' - "unrestricted permutation"
sens_attribute = 13  # column index of real dataset attribute to permute (1-13 for cheese dataset)


# for simulations (random data generation)
simulation_mode = 'random_equal'  # how to simulate 
                                  # options: 'random_equal' - equal sigma_i between assessors
                                  #          'random_unequal' - unequal sigma_i between assessors, sampled from uniform distribution
n_ass = 12  # I
n_samp = 14  # J
n_rep = 2   # K

#############################


# read command line arguments. options specified in command line will overwrite options hard-coded above
spec = matrix(c(
  'nsim', 'n', 2, "double",    # number of iterations
  'nass', 'i', 2, "double",    # number of assessors
  'nsam', 'j', 2, "double",    # number of samples
  'nrep', 'k', 2, "double",    # number of replicates
  'mode', 'm', 2, "character", # for permutations: use real data?
  'ratt', 'u', 2, "integer",   # column to use in real dataset 
  'ptyp', 't', 2, "character", # permutation type  ('ass', 'repl', 'assrepl', 'all')
  'simm', 's', 2, "character"   # simulation mode
), byrow=TRUE, ncol=4)
opt = getopt(spec)

aliases = c("nsim"="n_sim", "mode"="mode.", "nass"="n_ass", "nsam"="n_samp", "nrep"="n_rep",
            "ratt"='sens_attribute', "ptyp"="perm_type", "simm"="simulation_mode")

# overwrite options if set
for(arg_ in names(opt)[!names(opt) %in% c("ARGS", "taskid")])
  assign(aliases[names(aliases)==arg_], opt[[arg_]])






if(mode.=='permutation'){
  cheese = read.csv('cheese.csv')
  cheese_attributes = names(cheese)[-(1:3)]
  if(sens_attribute < 1 | sens_attribute > 13) stop ('attibute index out of range')
  
  sens_attribute = cheese_attributes[sens_attribute]
  real_dat = data.frame(intensity=cheese[,sens_attribute],
                        assessor=as.factor(cheese$Assessor),
                        sample=as.factor(cheese$Sample),
                        repl=as.factor(cheese$Replicates))

  n_ass = length(unique(real_dat$assessor))
  n_samp = length(unique(real_dat$samp))
  n_rep = length(unique(real_dat$repl))
  perm = T 
  simulation_mode = ''
} else if (mode.=='simulation'){
  perm = F
  perm_type = ''
} else {
  stop('set mode to simulation or permutation')
}

cat('###   settings   ###\n')
if(mode.=='permutation') {
  cat(paste0('permuting cheese data, attribute: ', sens_attribute, 
             '\npermutation mode: ', perm_type))
  
} else { 
  cat(paste0('simulating data', 
             '\nn_ass\t\t', n_ass, 
             '\nmode\t\t', simulation_mode, 
             '\nn_samp\t\t', n_samp, 
             '\nn_rep\t\t', n_rep))
}
cat(paste0('\nnumber of iterations: ', n_sim))  
print(  perm)
print(mode.)

outfile = paste0('out_', mode., '_', 
                 ifelse(mode.=='permutation', paste0(perm_type, '_cheese-', sens_attribute), 
                        paste0(simulation_mode, '_I', n_ass, '_J', n_samp, '_K', n_rep)),
                 '_', starttime, '.Rds')
cat(paste0('\nwriting to: ', paste0(getwd(), '/', outfile), '\n'))


# tables holding permutation / simulation results
f_table = p_table = matrix(nrow=n_sim, ncol=9, dimnames=list(1:n_sim, 
      c('raw', 'std', 'tenb', 'add', 'mult', 'MAM.AOV_default', 'MAM_adj', 'MAM.AOV_nonadj', 'MAM_nonadj')))

# additional outputs of the MAM
mam_extra = matrix(nrow=n_sim, ncol=5, dimnames=list(NULL, c('p_scaling_default', 'p_scaling_nonadj', 'p_disagree_default', 'p_disagree_nonadj', 'n_negatives')))


# permutation: use real (cheese) data
if(perm)
    sim_dat = real_dat  # inject real data


for(ns in 1:n_sim){  # simulation / permutation iterations
  if (ns %% 10 == 0) cat('iteration ', ns, '\n')
  status = F
  
  while(!status){
    
    if(perm & ns > 1) {  # real data: perform permutation
                         # 1st row of results table is UNpermuted data!
      if (perm_type == 'ass'){  
        # permute within assessors ("restricted permutation")
        # underlying model: y_{ijk} = mu + a_i + epsilon_{ijk}; epsilon ~ N(0,sigma_i^2)
        for (ass in unique(sim_dat$assessor)){
          ass_idx = which(sim_dat$assessor == ass)
          ass_idx_perm = sample(ass_idx, length(ass_idx), replace=F)
          sim_dat[ass_idx, 'intensity'] = sim_dat[ass_idx_perm, 'intensity']
        }
      } else if (perm_type == 'assrepl'){  
        # keep replicates and assessors together [not examined in the paper; not a feasible null hypothesis]
        # $y_{ijk} = mu + a_i + g_{ij} + epsilon_{ijk}; epsilon_{ijk} ~ N(0,\sigma_{ij}^2
        for (ass in unique(sim_dat$assessor)){
          ass_idx = which(sim_dat$assessor == ass)
          idx_ = sample(n_samp)
          ass_idx_perm = c(rbind(2*idx_-1, 2*idx_))
          sim_dat[ass_idx, 'intensity'] = sim_dat[ass_idx[ass_idx_perm], 'intensity']   # permute but keep replicates and assessors together
        } 
      } else if (perm_type == 'all'){  
        # permutation across all values (completely random, "unrestricted permutation")
        # underlying model: y_{ijk} = mu + epsilon_{ijk}; epsilon ~ N(0,sigma^2)
        sim_dat$intensity = sample(sim_dat$intensity) 
      } else stop('invalid permutation type')
    }
    
    if(mode.=='simulation') {  
      if(simulation_mode=='random_equal') {  # mu = 5, sigma=1.5
        sim_dat = data.frame(intensity=rnorm(n_ass*n_samp*n_rep, mean=5, sd=1.5),
                             assessor=as.factor(rep(1:n_ass, each=n_samp*n_rep)),
                             sample=as.factor(rep(1:n_samp, each=n_rep)),
                             repl=1:n_rep)
      
      } else if(simulation_mode=='random_unequal') {  # mu = 5, sigma=unif(0.1,1)
        sim_dat = data.frame(intensity=unlist(sapply(runif(n_ass, 0.1, 1), function(x) rnorm(n_samp*n_rep, mean=5, sd=x), simplify=F)),
                             assessor=as.factor(rep(1:n_ass, each=n_samp*n_rep)),
                             sample=as.factor(rep(1:n_samp, each=n_rep)),
                             repl=1:n_rep)
      } else {stop('invalid simulation mode')}
    }
    
    
    
    out_aov = scaling_fit_aov(sim_dat)
    if(!out_aov$status) next
    status = T
    out_mam = fit_MAM(sim_dat)
    f_table[ns,] = c(out_aov$aov_f, out_mam$mam_f)
    p_table[ns,] = c(out_aov$aov_p, out_mam$mam_p)
    mam_extra[ns,] = out_mam$mam_extra
    
  } # end while
  

  if ((10*ns / n_sim) %% 1 == 0 | (10*ns / (n_sim-1)) %% 1 == 0 | ns == n_sim){  # save at every 10% progress
    cat(paste('checkpoint at', ns, 'simulations. writing to disk:', outfile), '\n')
    sim_results = list(f_table=f_table, p_table=p_table, mam_extra=mam_extra,
                       parameters=list(mode.=mode., perm=perm, perm_type=perm_type, 
                                       sens_attribute=sens_attribute,
                                       simulation_mode=simulation_mode, n_ass=n_ass, n_samp=n_samp, n_rep=n_rep, n_sim=n_sim))
    saveRDS(sim_results, paste0(outfile)) 
  }
}  


print(paste('endtime:', format(Sys.time(), "%Y-%m-%d_%H.%M.%S")))
