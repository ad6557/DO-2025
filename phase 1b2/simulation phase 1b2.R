library(tidyverse)
library(survival)
library(Mediana)
library(ldbounds)
library(nph)
library(reshape2)
library(ggplot2)
library(parallel)
library(doParallel)
library(dplyr)
library(RBesT)
library(dplyr)


source("functions.R")


# comparable with seamless phase 2/3
p_ORR_h = 0.10
p_TOX_h = 0.30
p_ORR_l = 0.10
p_TOX_l = 0.30
p_ORR_c = 0.10
p_TOX_c = 0.10

N1 <- 50
N2 <- 370
N = N1+N2
n_arm_stg1 <- 3
n_arm_stg2 <- 2
##calculate nominal p-values using alpha spending function ----for OS superiority
alpha0 =0.025
times.OS=c(0.5,1);T = length(times.OS)
sigcut =rep(0,T); pcut=sigcut
sigcut=ldBounds(t=times.OS, alpha=alpha0, iuse=c(1), sides = 1)$upper.bounds ## 10/19/2024 Haiyang Sheng
#iuse = 1 OR 2 : O'Brien Fleming OR Pocock type boundaries respectively. 
##critical values; nominal values
pcut.OS=1-as.numeric(lapply(sigcut, pnorm))

# survival data
lambda_c = 0.06
HR = c(1,1) # c(HR_low,HR_high)

# Desired latent correlations:
# tox, eff, survival
rho_tox_eff <- 0.6
rho_tox_surv <- -0.5
rho_eff_surv <- -0.5
# Construct correlation matrix (must be symmetric and positive definite)
Sigma <- matrix(c(1, rho_tox_eff, rho_tox_surv,
                  rho_tox_eff, 1, rho_eff_surv,
                  rho_tox_surv, rho_eff_surv, 1), nrow = 3)

# closed form phi
phi_h = gen_phi(p_ORR_h,p_TOX_h,rho_tox_eff)
phi_l = gen_phi(p_ORR_l,p_TOX_l,rho_tox_eff)
phi_c = gen_phi(p_ORR_c,p_TOX_c,rho_tox_eff)


######################## step 1 boundary ##########################
set.seed(2025)
p_ORR_0 = 0.10 # meaningful efficacy under H0
p_TOX_0 = 0.30 # tolerable toxicity under H0
p_ORR_1 = 0.25 # meaningful efficacy under Ha
p_TOX_1 = 0.10 # tolerable toxicity under Ha
# Stopping boundaries for high and low dose
phi_0=gen_phi(p_ORR_0,p_TOX_0,rho_tox_eff)
phi_1=gen_phi(p_ORR_1,p_TOX_1,rho_tox_eff)
source("BOP2.R")
#stopBoundary_mat_opt

####################### step 2 boundary ###########################
eta_ORR = 0.05 ## Superiority margin of ORR
eta_TOX = 0.05 ## Superiority margin of TOX
## Selection error
alpha_L = 0.3 ## Not Recommend Low  when Low is better
alpha_H = 0.3## Not Recommend High when High is better

source("pick the winner.R")
#theta2_opt
#ORR_boundary_interim_tbl_list
#TOX_boundary_interim_tbl_list

####################### simulation starts #########################

# PoC dual criteria
decision = decision2S(0.8, 0, lower.tail=FALSE) #Pr(dORR>0)>80%
decision.1 = decision2S(0.5, 0.1, lower.tail=FALSE) #Pr(dORR>0.1)>50%

M=5000;#number of simulations

cols= c("N","DODII IA1","DODII IA2","PoC","phase III","IA decision","decision",
        "DS time","IA time","FA time")
output_list <- list()

fup = 1; # minimum follow up time for subject at dose selection analysis in Stage 1 (month)
D=2; # additional 2 more month for DBL and analysis
E.OS = 390 ; #total number of events in phase 3

for (HH in c("H1")){ # ###type I error: H0 ; power : H1
  
  output= as.data.frame(matrix(0, M,length(cols)));colnames(output)=cols
  output$N = N
  
  if (HH == "H1"){
    p_ORR = c(p_ORR_c,p_ORR_l,p_ORR_h)
    p_TOX = c(p_TOX_c,p_TOX_l,p_TOX_h)
    lambda = c(lambda_c, HR[1]*lambda_c, HR[2]*lambda_c)
  } else if (HH == "H0"){
    p_ORR = c(p_ORR_c,p_ORR_c,p_ORR_c)
    p_TOX = c(p_TOX_c,p_TOX_c,p_TOX_c)
    lambda = c(lambda_c, lambda_c, lambda_c)
  }
  
  for (Iter in 1:M){
    print(paste0("Iter: ",Iter))
    seed= Iter*10 -1 ##set seed to get reproducible results
    select = 0
    
    ### simulate ORR and TOX and OS data for phase 2 and phase 3
    
    data = sim_data(N1,N2,lambda,seed,p_ORR,p_TOX,Sigma) ##simulated data for OS
    data$flup1 = NA; #follow up time at IA1 (dose selection analysis)
    #fup is the minimum follow up time at dose selection DBL;
    #assume enroll 15 subjects in total (aka 5 subjects per arm) per month
    mfup =fup +ceiling((N1*n_arm_stg1)/15)-1 #maximum follow up at dose selection
    data$flup1[which(data$Group==0& data$Stage==1)]= sample(rep(fup:mfup,each=5),N1) # enrollment 10 people per month per arm
    data$flup1[which(data$Group==1& data$Stage==1)]= sample(rep(fup:mfup,each=5),N1)
    data$flup1[which(data$Group==2& data$Stage==1)]= sample(rep(fup:mfup,each=5),N1)
    data$dif = mfup-data$flup1
    
    # IA1 at DODII
    data_S1 <- data %>%
      filter(Stage == 1) %>%
      group_by(Group) %>%
      arrange(dif) %>%
      slice_head(n = 25) %>% # ties are kept in the order they appear in the data
      ungroup()
    n_High_ORR = sum(subset(data_S1, Group == 2)$orr)
    n_Low_ORR = sum(subset(data_S1, Group == 1)$orr)
    n_High_TOX = sum(subset(data_S1, Group == 2)$tox)
    n_Low_TOX = sum(subset(data_S1, Group == 1)$tox)
    BOP2_h = BOP2_gatekeeper(n_High_ORR,n_High_TOX,1,stopBoundary_mat_opt)
    BOP2_l = BOP2_gatekeeper(n_Low_ORR,n_Low_TOX,1,stopBoundary_mat_opt)
    if(BOP2_h==0&BOP2_l==0){
      dose_bop1 = 0
      output[Iter, "DS time"] = max(data_S1$dif)+fup+D # 7
      next
    } else if (BOP2_h==0&BOP2_l==1){
      dose_bop1 = 1
      output[Iter, "DS time"] <- mfup+D # 12
    } else if (BOP2_h==1&BOP2_l==0){
      dose_bop1 = 2
      output[Iter, "DS time"] <- mfup+D # 12
    } else if (BOP2_h==1&BOP2_l==1){
      dose_bop1 = pick_winner(n_High_ORR, n_Low_ORR, n_High_TOX, n_Low_TOX, 1, ORR_boundary_interim_tbl_list, TOX_boundary_interim_tbl_list)
      output[Iter, "DS time"] <- mfup+D # 12
    }
    output[Iter, "DODII IA1"] = sum(dose_bop1) #1:low; 2:high; 3:both
    
    # IA2 at DODII
    data_S2 = subset(data, Stage == 1)
    n_High_ORR2 = if (2 %in% dose_bop1) sum(subset(data_S2, Group == 2)$orr) else NA
    n_Low_ORR2  = if (1 %in% dose_bop1) sum(subset(data_S2, Group == 1)$orr) else NA
    n_High_TOX2 = if (2 %in% dose_bop1) sum(subset(data_S2, Group == 2)$tox) else NA
    n_Low_TOX2  = if (1 %in% dose_bop1) sum(subset(data_S2, Group == 1)$tox) else NA
    if(is.na(n_High_ORR2)){
      BOP2_h2 = 0
      BOP2_l2 = BOP2_gatekeeper(n_Low_ORR2,n_Low_TOX2,2,stopBoundary_mat_opt)
    }else if(is.na(n_Low_ORR2)){
      BOP2_h2 = BOP2_gatekeeper(n_High_ORR2,n_High_TOX2,2,stopBoundary_mat_opt)
      BOP2_l2 = 0
    }else{
      BOP2_h2 = BOP2_gatekeeper(n_High_ORR2,n_High_TOX2,2,stopBoundary_mat_opt)
      BOP2_l2 = BOP2_gatekeeper(n_Low_ORR2,n_Low_TOX2,2,stopBoundary_mat_opt)
    }
    
    # Final decision logic
    if (BOP2_h2 == 0 & BOP2_l2 == 0) {
      dose_bop2 <- 0
      next
    } else if (BOP2_h2 == 0 & BOP2_l2 == 1) {
      dose_bop2 <- 1
    } else if (BOP2_h2 == 1 & BOP2_l2 == 0) {
      dose_bop2 <- 2
    } else if (BOP2_h2 == 1 & BOP2_l2 == 1) {
      dose_bop2 <- pick_winner(n_High_ORR2, n_Low_ORR2, n_High_TOX2, n_Low_TOX2,
                               2, ORR_boundary_interim_tbl_list, TOX_boundary_interim_tbl_list)
    }
    
    if(length(dose_bop2)==2 & n_High_ORR2>n_Low_ORR2){
      dose_bop2 <- 2
    }else if (length(dose_bop2)==2 & n_High_ORR2<=n_Low_ORR2){
      dose_bop2 <- 1
    }
    
    select = dose_bop2
    output[Iter, "DODII IA2"] = select
    if(select==0){next}
    
    ### PoC
    orr_low <- mean(subset(data_S2, Group == select)$orr) # ORR for treatment arm
    orr_c <- mean(subset(data_S2, Group == 0)$orr) # ORR for control
    
    nf.prior <- mixbeta(nf.prior = c(1,0.5,0.5)) # non-informative prior
    
    ss_trt <- N1
    ss_control <- N1
    
    post_c <- postmix(nf.prior, n = ss_control, r = sum(subset(data_S2, Group == 0)$orr)) #control: non-informative
    post_t <- postmix(nf.prior, n = ss_trt, r = sum(subset(data_S2, Group == select)$orr)) #trt DL2: non-informative
    p.t <- decision(post_t, post_c)+decision.1(post_t, post_c) #non-informative post-t
    
    output[Iter, "PoC"] = p.t
    if(p.t==0){next}
    
    
    ### phase 3
    data_s <- subset(data, (Group == 0 | Group == select) & Stage == 2, select = -c(orr, tox))
    data_s$Group[which(data_s$Group==select)]=1
    
    
    data_s$dif=NA # the difference between subject enroll time from the first enrollment (t=0)
    #---so firstly 10 per arm until dose selection (for fup month); 
    #then 1 month followup 2 months read out
    #then enroll 20 subject per arm until reach N subjects per arm
    remd = N2%%20; quot= N2%/%20#remainder and the quotient :
    data_s$dif[which(data_s$Group==0& data_s$Stage==2)]= sample(c(12+rep(1:quot-1,each=20),+rep(12+quot,remd)),N2)
    data_s$dif[which(data_s$Group==1& data_s$Stage==2)]= sample(c(12+rep(1:quot-1,each=20),+rep(12+quot,remd)),N2)
    table(data_s$dif)
    data_s$ref.OS = data_s$Time+ data_s$dif #calender time from the first enrollment
    sort = data_s[order(data_s$ref.OS),];
    sort$look.OS = 0;
    tcut1= sort$ref.OS[ceiling(E.OS*0.5)];
    output[Iter,c("IA time")] = tcut1
    sort$look.OS[which(sort$ref.OS <= tcut1)] =1; 
    #final analysis of OS when E.OS occurs
    tcut2= sort$ref.OS[ceiling(E.OS)];
    output[Iter,c("FA time")] = tcut2
    sort$look.OS[which(sort$look.OS ==0 & sort$ref.OS <= tcut2)] =2;
    
    #############################################################################################
    ###################
    ##at IA OS look (50% OS events)
    dat1 = sort[c("Time","look.OS","Group","dif")]; 
    colnames(dat1) =c("Time","look","Group","dif")
    dat1$ind= ifelse(dat1$look %in% 1 , 1, 0 ); sum(dat1$ind)
    dat1$flup=tcut1-dat1$dif; 
    dat1$Time[which(dat1$ind==0)] = dat1$flup[which(dat1$ind==0)];
    dat1 <- dat1[dat1$Time >= 0, ]
    pval1.OS = surv.test(dat1, altern = "greater") ## Superiority
    
    #############################################################################################
    ###################
    ##at OS final look (100% OS events)
    dat2 = sort[c("Time","look.OS","Group","dif")]; 
    colnames(dat2) =c("Time","look","Group", "dif")
    dat2$ind= ifelse(dat2$look %in% 1:2 , 1, 0 ); sum(dat2$ind)
    dat2$flup=tcut2-dat2$dif; head(dat2);##add censored values for censored case
    dat2$Time[which(dat2$ind==0)] = dat2$flup[which(dat2$ind==0)];
    pval2.OS = surv.test(dat2, altern = "greater") ## Superiority
    
    ##SAVE
    
    Allps_noadj= c(pval1.OS,pval2.OS)
    output[Iter,c("IA decision")]= as.numeric(Allps_noadj[1]<pcut.OS[1])
    p3 =as.numeric(sum(Allps_noadj<pcut.OS)>0)
    output[Iter,c("phase III")] = p3
    output[Iter,c("decision")]= p3 * (select > 0) * (p.t > 0) # suppose to be the sampe as p3
  } #end of M
  output_list[[HH]] <- output
} #end of HH


# overall power
mean(output_list[["H1"]][["decision"]])
# overall type I error
mean(output_list[["H0"]][["decision"]])

# phase 3 power
sum(output_list[["H1"]][["phase III"]])/sum(output_list[["H1"]][["PoC"]] != 0)
# phase 3 type I error
sum(output_list[["H0"]][["phase III"]])/sum(output_list[["H0"]][["PoC"]] != 0)




output_suc <- output %>%
  filter(decision != 0) %>%
  mutate(
    ess = case_when(
      `IA decision` == 0 ~ 890,
      `IA decision` == 1 & (`IA time` + 2) >= 30 ~ 890,
      `IA decision` == 1 & (`IA time` + 2) >= 29 ~ 870,
      `IA decision` == 1 & (`IA time` + 2) >= 28 ~ 830,
      `IA decision` == 1 & (`IA time` + 2) >= 27 ~ 790,
      `IA decision` == 1 & (`IA time` + 2) >= 26 ~ 750,
      `IA decision` == 1 & (`IA time` + 2) >= 25 ~ 710,
      `IA decision` == 1 & (`IA time` + 2) >= 24 ~ 670,
      TRUE ~ NA_real_
    )
  )
ess1 = mean(output_suc$ess)
duration = sum(output_suc$`FA time`*(1-output_suc$`IA decision`)+output_suc$`IA time`*output_suc$`IA decision`)/dim(output_suc)[1]
duration1 = duration+2 # read out 2 months
duration1
mean(output_suc$`FA time`) # then 2 months read out follow it
mean(output_suc$`IA time`) # then 2 months read out follow it



output_fail = output[output$`decision`==0,]
output_fail$ess[output_fail$`DS time`==7] = 7*5*3 
output_fail$ess[output_fail$`DS time`==12&output_fail$`FA time` == 0] = 50*3 
output_fail$ess[output_fail$`FA time` != 0] = 890 
ess0 = mean(output_fail$ess)
idx <- which(output_fail$`FA time` == 0)
output_fail$`FA time`[idx] <- output_fail$`DS time`[idx]
output_fail$`FA time`[-idx] <- output_fail$`FA time`[-idx]+2 # read out 2 months
duration0 = mean(output_fail$`FA time`)
duration0


table(output_list[["H1"]][["DODII IA1"]])
table(output_list[["H0"]][["DODII IA1"]])

table(output_list[["H1"]][["DODII IA2"]])
table(output_list[["H0"]][["DODII IA2"]])

table(output_list[["H1"]][["PoC"]])
table(output_list[["H0"]][["PoC"]])

table((output[["PoC"]]>0))/M
table(output[["DODII IA2"]],(output[["PoC"]]>0))/M


# when HR_low = 1, HR_high = 0.7
# overall power
mean(output[["decision"]]&output[["DODII IA2"]]==2)
# phase 3 power
sum(output[["decision"]]==1&output[["DODII IA2"]]==2)/sum(output[["PoC"]] != 0)





HR_string <- paste(HR, collapse = "_")
save.image(file = paste0("S02_N2_",N2,"_HR_",HR_string,".RData"))