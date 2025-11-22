# I didn't write as a function, as %dopar% is hard to deal with.

## External parameters
N =N1 # maximum sample size
L =2 # The number of stages
n = N/L # The number of new patients enrolled at each stage
# p_ORR_0 = 0.3 # meaningful efficacy under H0
# p_TOX_0 = 0.75 # tolerable toxicity under H0
# p_ORR_1 = 0.45 # meaningful efficacy under Ha
# p_TOX_1 = 0.60 # tolerable toxicity under Ha
n.sim =10000 # number of simulation
# phi: Strength of association between ORR and TOX
######### Pick up BOPII stopping Boundary (Grid search tunning parameters)
grids.num_lmd = 1000; # grid resolution
grids.num_nu=1000;# grid resolution
lmd_L = 0.1; # Lower bound for lambda
lmd_H = 0.9; # Upper bound for lambda
nu_L = 0.1; # Lower bound for nu
nu_H=2; # Upper bound for nu
# To get the marginal of ORR and TOX in a multinomial vector
collapse_ORR = c(1,1,0,0) ## collapse pi1, pi2
collapse_TOX = c(1,0,1,0) ## collapse pi1, pi3
alpha = 0.05 # Type I error control;

num_core = detectCores()-1

## Construct the grid
lmd_vec = seq(lmd_L,lmd_H,length.out=grids.num_lmd)
nu_vec = seq(nu_L,nu_H,length.out = grids.num_nu)

## Simulation for L stages
pi_with_phi_H0 = gen_pis(p_ORR_0,p_TOX_0,phi=phi_0)
pi_with_phi_H1 = gen_pis(p_ORR_1,p_TOX_1,phi=phi_1)
sim_array_H0 = Trial_sim(pi_with_phi_H0,n.sim=n.sim,L,n.interim=n)
sim_array_H1 = Trial_sim(pi_with_phi_H1,n.sim=n.sim,L,n.interim=n)
cl=makeCluster(num_core)
registerDoParallel(cl)

lambda_nu_opt = as.matrix(expand.grid(lmd_vec,nu_vec)) # (1000*1000) x 2
## Stopping boundaries for 
lambda_nu_opt_list=split(lambda_nu_opt,row(lambda_nu_opt))
stopBoundary_mat_list= parLapply(cl,lambda_nu_opt_list,
                                 stopBoundary_cal,L=L,n=n, p_ORR_0=p_ORR_0,p_TOX_0=p_TOX_0,
                                 collapse_ORR=collapse_ORR,collapse_TOX=collapse_TOX,
                                 pi_with_phi = pi_with_phi_H0,n0=1)

stopBoundary_mat_list_unique_raw = unique(stopBoundary_mat_list)
stopBoundary_mat_list_unique = lapply(stopBoundary_mat_list_unique_raw,function(x){
  if(!anyNA(x)){
    return(x)
  }
})
stopBoundary_mat_list_unique = stopBoundary_mat_list_unique[lapply(stopBoundary_mat_list_unique,length)>0]


t1e_vec = foreach(list_i=1:length(stopBoundary_mat_list_unique),.combine = c,.export = c("power_cal")) %dopar%
  mean(apply(sim_array_H0,MARGIN = 3,power_cal,stopBoundary_mat=stopBoundary_mat_list_unique[[list_i]]))
power_vec = foreach(list_i=1:length(stopBoundary_mat_list_unique),.combine = c,.export = c("power_cal")) %dopar%
  mean(apply(sim_array_H1,MARGIN = 3,power_cal,stopBoundary_mat=stopBoundary_mat_list_unique[[list_i]]))
param.opt = which(power_vec==max(power_vec[t1e_vec<=alpha]))
stopBoundary_mat_opt = stopBoundary_mat_list_unique[[param.opt]]