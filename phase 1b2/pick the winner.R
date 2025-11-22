phi_High = vector("numeric", length = 8)
phi_Low = vector("numeric", length = 8)
# The ORR and TOX assumptions for those 8 scenarios 
## H_1
p_ORR_High_H1 = p_ORR_1+eta_ORR;p_TOX_High_H1 = p_TOX_1
p_ORR_Low_H1 = p_ORR_1;p_TOX_Low_H1 = p_TOX_1
phi_High[1] = gen_phi(p_ORR_High_H1,p_TOX_High_H1,rho_tox_eff)
phi_Low[1] = gen_phi(p_ORR_Low_H1,p_TOX_Low_H1,rho_tox_eff)

## H_2
p_ORR_High_H2 = p_ORR_1+eta_ORR;p_TOX_High_H2 = p_TOX_1
p_ORR_Low_H2 = p_ORR_1;p_TOX_Low_H2 = p_TOX_1-eta_TOX
phi_High[2] = gen_phi(p_ORR_High_H2,p_TOX_High_H2,rho_tox_eff)
phi_Low[2] = gen_phi(p_ORR_Low_H2,p_TOX_Low_H2,rho_tox_eff)

## H_3
p_ORR_High_H3 = p_ORR_1;p_TOX_High_H3 = p_TOX_1
p_ORR_Low_H3 = p_ORR_1;p_TOX_Low_H3 = p_TOX_1
phi_High[3] = gen_phi(p_ORR_High_H3,p_TOX_High_H3,rho_tox_eff)
phi_Low[3] = gen_phi(p_ORR_Low_H3,p_TOX_Low_H3,rho_tox_eff)

## H_4
p_ORR_High_H4 = p_ORR_1;p_TOX_High_H4 = p_TOX_1
p_ORR_Low_H4 = p_ORR_1;p_TOX_Low_H4 = p_TOX_1-eta_TOX
phi_High[4] = gen_phi(p_ORR_High_H4,p_TOX_High_H4,rho_tox_eff)
phi_Low[4] = gen_phi(p_ORR_Low_H4,p_TOX_Low_H4,rho_tox_eff)

## H_5
p_ORR_High_H5 = p_ORR_0;p_TOX_High_H5 = p_TOX_0
p_ORR_Low_H5 = p_ORR_0;p_TOX_Low_H5 = p_TOX_0
phi_High[5] = gen_phi(p_ORR_High_H5,p_TOX_High_H5,rho_tox_eff)
phi_Low[5] = gen_phi(p_ORR_Low_H5,p_TOX_Low_H5,rho_tox_eff)

## H_6
p_ORR_High_H6 = p_ORR_0;p_TOX_High_H6 = p_TOX_0
p_ORR_Low_H6 = p_ORR_0;p_TOX_Low_H6 = p_TOX_1
phi_High[6] = gen_phi(p_ORR_High_H6,p_TOX_High_H6,rho_tox_eff)
phi_Low[6] = gen_phi(p_ORR_Low_H6,p_TOX_Low_H6,rho_tox_eff)

## H_7
p_ORR_High_H7 = p_ORR_1;p_TOX_High_H7 = p_TOX_0
p_ORR_Low_H7 = p_ORR_0;p_TOX_Low_H7 = p_TOX_0
phi_High[7] = gen_phi(p_ORR_High_H7,p_TOX_High_H7,rho_tox_eff)
phi_Low[7] = gen_phi(p_ORR_Low_H7,p_TOX_Low_H7,rho_tox_eff)

## H_8
p_ORR_High_H8 = p_ORR_1;p_TOX_High_H8 = p_TOX_0
p_ORR_Low_H8 = p_ORR_0;p_TOX_Low_H8 = p_TOX_1
phi_High[8] = gen_phi(p_ORR_High_H8,p_TOX_High_H8,rho_tox_eff)
phi_Low[8] = gen_phi(p_ORR_Low_H8,p_TOX_Low_H8,rho_tox_eff)


### Simulation for Comparison
for(case_i in 1:8){
  p_ORR_High_temp = get(paste0('p_ORR_High_H',case_i))
  p_TOX_High_temp = get(paste0('p_TOX_High_H',case_i))
  p_ORR_Low_temp = get(paste0('p_ORR_Low_H',case_i))
  p_TOX_Low_temp = get(paste0('p_TOX_Low_H',case_i))
  pi_with_phi_High_temp = gen_pis(p_ORR_High_temp,p_TOX_High_temp,phi=phi_High[case_i])
  pi_with_phi_Low_temp = gen_pis(p_ORR_Low_temp,p_TOX_Low_temp,phi=phi_Low[case_i])
  sim_array_Low_temp = Trial_sim(pi_with_phi_Low_temp,n.sim=n.sim,L,n.interim=n)
  sim_array_High_temp = Trial_sim(pi_with_phi_High_temp,n.sim=n.sim,L,n.interim=n)
  sim_array_High_Low_temp = array(NA,dim = c(4,L,n.sim,2)) #[4,2,10000,2]
  sim_array_High_Low_temp[,,,1] = sim_array_High_temp
  sim_array_High_Low_temp[,,,2] = sim_array_Low_temp
  assign(paste0("sim_array_High_Low_Compare_H",case_i),sim_array_High_Low_temp)
  assign(paste0('pi_with_phi_High_H',case_i),pi_with_phi_High_temp) # [pi_d1,pi_d2,pi_d3,pi_d4]
  assign(paste0('pi_with_phi_Low_H',case_i),pi_with_phi_Low_temp) # [pi_d1,pi_d2,pi_d3,pi_d4]
}

## List of 8 scienarios along with their prior assumptions
sim_array_High_Low_Compare_list = list(sim_array_High_Low_Compare_H1,
                                       sim_array_High_Low_Compare_H2,
                                       sim_array_High_Low_Compare_H3,
                                       sim_array_High_Low_Compare_H4,
                                       sim_array_High_Low_Compare_H5,
                                       sim_array_High_Low_Compare_H6,
                                       sim_array_High_Low_Compare_H7,
                                       sim_array_High_Low_Compare_H8)

pi_with_phi_Low_list = list(pi_with_phi_Low_H1,
                            pi_with_phi_Low_H2,
                            pi_with_phi_Low_H3,
                            pi_with_phi_Low_H4,
                            pi_with_phi_Low_H5,
                            pi_with_phi_Low_H6,
                            pi_with_phi_Low_H7,
                            pi_with_phi_Low_H8)

pi_with_phi_High_list = list(pi_with_phi_High_H1,
                             pi_with_phi_High_H2,
                             pi_with_phi_High_H3,
                             pi_with_phi_High_H4,
                             pi_with_phi_High_H5,
                             pi_with_phi_High_H6,
                             pi_with_phi_High_H7,
                             pi_with_phi_High_H8)







## Calculating the OC for different theta2.
theta2_vec = seq(0.51,0.99,by=0.01)
system.time({
  theta2_grid_raw=foreach(theta2 = theta2_vec,.combine = rbind) %:%
    foreach(case_id = 1:4, .combine = rbind) %dopar%
    stopBoundary_cal_compare(theta2_ORR=rep(theta2,L),theta2_TOX = rep(theta2,L),case_id=case_id,
                             sim_array_High_Low_Compare_list=sim_array_High_Low_Compare_list,
                             pi_with_phi_High_list=pi_with_phi_High_list,
                             pi_with_phi_Low_list=pi_with_phi_Low_list,stopBoundary_mat_opt = stopBoundary_mat_opt,n=n)
})


# Tabulating
theta2_tbl=tibble(
  case = theta2_grid_raw[,10],
  theta2 = as.numeric(theta2_grid_raw[,8]),
  Recommend_High_only = as.numeric(theta2_grid_raw[,1]),
  Recommend_Low_only = as.numeric(theta2_grid_raw[,2]),
  Recommend_Both = as.numeric(theta2_grid_raw[,3]),
  Recommend_None = as.numeric(theta2_grid_raw[,4]),
  n_sample_size = as.numeric(theta2_grid_raw[,5]),
  n_sample_size_High_dose = as.numeric(theta2_grid_raw[,6]),
  n_sample_size_Low_dose = as.numeric(theta2_grid_raw[,7])
)


#### Optimization criteria to select the optimal theta2
theta2_tbl_alpha_L_controlled=
  theta2_tbl%>%filter(case==4&Recommend_High_only+Recommend_None<=alpha_L)%>%
  dplyr::select(theta2)

theta2_tbl_alpha_H_controlled=theta2_tbl%>%filter(case==1&Recommend_Low_only+Recommend_None<=alpha_H)%>%
  dplyr::select(theta2)


theta2_cands_tbl= theta2_tbl_alpha_H_controlled%>%inner_join(theta2_tbl_alpha_L_controlled,
                                                             by=c("theta2"))%>%
  distinct()


theta2_powers_tbl = theta2_cands_tbl%>%left_join(theta2_tbl,by=c("theta2"))%>%
  filter(case%in%c(1,4))


beta_H_tibble = theta2_powers_tbl%>%
  filter(case==1)%>%
  dplyr::select(theta2,Recommend_High_only)

beta_L_tibble = theta2_powers_tbl%>%
  filter(case==4)%>%
  dplyr::select(theta2,Recommend_Low_only)


beta_H_tibble$Recommend_Low_only = beta_L_tibble$Recommend_Low_only

beta_HL_tibble = beta_H_tibble%>%
  group_by(theta2)%>%
  mutate(min_power = min(Recommend_High_only,Recommend_Low_only))


theta2_opt = beta_HL_tibble%>%
  filter(min_power==max(beta_HL_tibble$min_power))%>%
  dplyr::select(theta2)%>%tail(1)%>%pull()







theta2_ORR = rep(theta2_opt,L); 
theta2_TOX = rep(theta2_opt,L)

##############################################################################
collapse_ORR = c(1,1,0,0) ## collapse pi1, pi2
collapse_TOX = c(1,0,1,0) ## collapse pi1, pi3

for(case_i in 1:4){
  p_ORR_High_temp = get(paste0('p_ORR_High_H',case_i))
  p_TOX_High_temp = get(paste0('p_TOX_High_H',case_i))
  p_ORR_Low_temp = get(paste0('p_ORR_Low_H',case_i))
  p_TOX_Low_temp = get(paste0('p_TOX_Low_H',case_i))
  pi_with_phi_High_temp = gen_pis(p_ORR_High_temp,p_TOX_High_temp,phi=phi_High[case_i])
  pi_with_phi_Low_temp = gen_pis(p_ORR_Low_temp,p_TOX_Low_temp,phi=phi_Low[case_i])
  assign(paste0('pi_with_phi_High_H',case_i),pi_with_phi_High_temp)
  assign(paste0('pi_with_phi_Low_H',case_i),pi_with_phi_Low_temp)
}




cl = makeCluster(num_core)
registerDoParallel(cl)


### Get integer boundary
ORR_boundary_interim_tbl_list = list()
TOX_boundary_interim_tbl_list = list()
for(L_i in 1:L){
  
  # pos_prob_compare_ORR_temp is P(ORR_h - ORR_l > eta_ORR | X)  
  pos_prob_ORR_combinations = foreach(X_H_i = (stopBoundary_mat_opt[L_i,2]+1):(L_i*n),.combine = cbind) %:%
    foreach(X_L_i=(stopBoundary_mat_opt[L_i,2]+1):(L_i*n),.combine=rbind) %dopar%
    compare_pos_cal(eta = eta_ORR,x_1 = c(X_H_i,0,0,L_i*n-X_H_i),x_2 = c(X_L_i,0,0,L_i*n-X_L_i),
                    collapse_type = collapse_ORR,pi_with_phi_1 = pi_with_phi_High_H3,
                    pi_with_phi_2 = pi_with_phi_Low_H3,n0 = 1)
  
  
  # pos_prob_compare_TOX_temp is P(TOX_h - TOX_l < eta_TOX | X)  
  pos_prob_TOX_combinations = 1-foreach(X_H_i = 0:(stopBoundary_mat_opt[L_i,3]-1),.combine = cbind) %:%
    foreach(X_L_i = 0:(stopBoundary_mat_opt[L_i,3]-1),.combine=rbind) %dopar%
    compare_pos_cal(eta = eta_TOX,x_1 = c(X_H_i,0,0,L_i*n-X_H_i),x_2 = c(X_L_i,0,0,L_i*n-X_L_i),
                    collapse_type = collapse_TOX,pi_with_phi_1 = pi_with_phi_High_H3,
                    pi_with_phi_2 = pi_with_phi_Low_H3,n0 = 1)
  
  
  
  
  pos_prob_ORR_combinations_tbl = tibble(
    num_High_ORR = rep((stopBoundary_mat_opt[L_i,2]+1):(L_i*n),each = length((stopBoundary_mat_opt[L_i,2]+1):(L_i*n))),
    num_LOW_ORR = rep((stopBoundary_mat_opt[L_i,2]+1):(L_i*n),length((stopBoundary_mat_opt[L_i,2]+1):(L_i*n))),
    ORR_decision = c(ifelse(pos_prob_ORR_combinations<1-theta2_ORR[L_i],"ORR_high_sig_not_better",
                            ifelse(pos_prob_ORR_combinations>theta2_ORR[L_i],"ORR_high_sig_better","Unclear")))
  )
  
  
  
  
  pic_ORR_temp = pos_prob_ORR_combinations_tbl%>%
    ggplot()+
    geom_rect(aes(xmin=num_LOW_ORR,
                  xmax=num_LOW_ORR+.75,
                  ymin=num_High_ORR,
                  ymax = num_High_ORR+.75,fill = ORR_decision))+
    scale_fill_manual(values =c("ORR_high_sig_better"="green",
                                "ORR_high_sig_not_better"="red",
                                "Unclear"="yellow"),
                      labels=c("High is better","High is not better","Unclear"),name="Decision for ORR")+
    ylab("Responses in High dose")+xlab("Responses in Low dose")+
    ggtitle(label = paste0("CutoffORR=",round(theta2_ORR[L_i],2)," CutoffTOX=",round(theta2_TOX[L_i],2)))
  
  
  
  pos_prob_TOX_combinations_tbl = tibble(
    num_High_TOX = rep(0:(stopBoundary_mat_opt[L_i,3]-1),each = length(0:(stopBoundary_mat_opt[L_i,3]-1))),
    num_LOW_TOX = rep(0:(stopBoundary_mat_opt[L_i,3]-1),length(0:(stopBoundary_mat_opt[L_i,3]-1))),
    TOX_decision = c(ifelse(pos_prob_TOX_combinations>theta2_TOX[L_i],"TOX_high_sig_not_worse",
                            ifelse(pos_prob_TOX_combinations<1-theta2_TOX[L_i],"TOX_high_sig_worse","Unclear")))
  )
  
  pic_TOX_temp =pos_prob_TOX_combinations_tbl%>%
    ggplot()+
    geom_rect(aes(xmin=num_LOW_TOX,
                  xmax=num_LOW_TOX+.75,
                  ymin=num_High_TOX,
                  ymax = num_High_TOX+.75,fill = TOX_decision))+
    scale_fill_manual(values =c("TOX_high_sig_not_worse"="green",
                                "TOX_high_sig_worse"="red",
                                "Unclear"="yellow"),
                      labels=c("High is not worse","High is worse","Unclear"),name="Decision for TOX")+
    ylab("Toxicities in High dose")+xlab("Toxicities in Low dose")+
    ggtitle(label = paste0("CutoffORR=",round(theta2_ORR[L_i],2)," CutoffTOX=",round(theta2_TOX[L_i],2)))
  filenames_ORR_temp = paste0("ORR_stopping_boundary_for_interim",L_i,"_same_theta_func_sample.pdf")
  filenames_TOX_temp = paste0("TOX_stopping_boundary_for_interim",L_i,"_same_theta_func_sample.pdf")
  ggsave(filename = filenames_ORR_temp,plot = pic_ORR_temp,width = 5,height = 5)
  ggsave(filename = filenames_TOX_temp,plot = pic_TOX_temp,width = 5,height = 5)
  ORR_boundary_interim_tbl_list[[L_i]] = pos_prob_ORR_combinations_tbl
  TOX_boundary_interim_tbl_list[[L_i]] = pos_prob_TOX_combinations_tbl
  
}

stopCluster(cl)