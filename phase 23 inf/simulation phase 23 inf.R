# dose expansion: monitor the efficacy of high dose and make a GNG decision
# - given p_ORR (no safety) in high dose from DODII
# if passes GNG
# seamless phase 2/3
# - generate TOX and ORR the same way in DODII using multinomial
# - generate correlated ORR/OS data

library(tidyverse)
library(survival)
library(Mediana)
library(ldbounds)
library(nph)
library(reshape2)
library(ggplot2)
library(dplyr)

source("functions.R")

## dose expansion
p_ORR_h = 0.10
p_TOX_h = 0.30
p_ORR_l = 0.10
p_TOX_l = 0.30
p_ORR_c = 0.10
p_TOX_c = 0.10

# Desired latent correlations:
# tox, eff, survival
rho_tox_eff <- 0.6
rho_tox_surv <- -0.5
rho_eff_surv <- -0.5
# Construct correlation matrix (must be symmetric and positive definite)
Sigma <- matrix(c(1, rho_tox_eff, rho_tox_surv,
                  rho_tox_eff, 1, rho_eff_surv,
                  rho_tox_surv, rho_eff_surv, 1), nrow = 3)
# survival data
lambda_c = 0.06
HR = c(1,1) # c(HR_low,HR_high)


## 3 arms in total for stage 1 and 2 arms in total for stage 2 
N1 <- 50
N2 <- 350
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





M=5000;#number of simulations

cols= c("N","expansion","select dose",
        "IA Dunnett","Dunnett","decision Dunnett","NoAdj","decision NoAdj","IA NoAdj",
        "DS time","IA time","FA time")
output_list <- list()

fup = 1; # minimum follow up time for subject at dose selection analysis in Stage 1 (month)
D=2; # additional 2 more month for DBL and analysis
E.OS = 390 ; #total number of events in both arms (Stage 1 + Stage 2)

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
    print(Iter)
    seed= Iter*10 -1 ##set seed to get reproducible results
    
    ### dose expansion
    # closed form phi
    phi_h = gen_phi(p_ORR[3],p_TOX[3],rho_tox_eff)
    pi_high <- gen_pis(p_ORR[3], p_TOX[3], phi_h)
    phi_l = gen_phi(p_ORR[2],p_TOX[2],rho_tox_eff)
    pi_low <- gen_pis(p_ORR[2], p_TOX[2], phi_l)
    # single arm: 20 people
    exp_h <- apply(rmultinom(20, 1, pi_high), 2, which.max)
    orr_h <- as.numeric(exp_h %in% c(1, 2))
    exp_l <- apply(rmultinom(20, 1, pi_low), 2, which.max)
    orr_l <- as.numeric(exp_l %in% c(1, 2))
    if(sum(orr_h)>2|sum(orr_l)>2){expanse=1}else{output[Iter,c("DS time")] = 6; next}
    output$expansion[Iter]=expanse
    
    ### simulate ORR and TOX and OS data for stage 1 and stage 2
    
    data = sim_data(N1,N2,lambda,seed,p_ORR,p_TOX,Sigma) ##simulated data for OS
    data$flup1 = NA; #follow up time at IA1 (dose selection analysis)
    #fup is the minimum follow up time at dose selection DBL;
    #assume enroll 30 subjects in total (aka 10 subjects per arm) per month
    mfup =fup +ceiling((N1*n_arm_stg1)/30)-1 #maximum follow up at dose selection
    output[Iter,c("DS time")] = 6+mfup+2 #dose expansion + maximum follow up + read out
    data$flup1[which(data$Group==0& data$Stage==1)]= sample(rep(fup:mfup,each=10),N1) # enrollment 10 people per month per arm
    data$flup1[which(data$Group==1& data$Stage==1)]= sample(rep(fup:mfup,each=10),N1)
    data$flup1[which(data$Group==2& data$Stage==1)]= sample(rep(fup:mfup,each=10),N1)
    
    ### Futility check
    data_S2 = subset(data, Stage == 1)
    # Determine if dose 1 or dose 2 should be dropped
    drop_dose1 <- (mean(data_S2[data_S2$Group == 1, "orr"]) <= 0.1 | 
                     mean(data_S2[data_S2$Group == 1, "tox"]) >= 0.3)
    drop_dose2 <- (mean(data_S2[data_S2$Group == 2, "orr"]) <= 0.1 | 
                     mean(data_S2[data_S2$Group == 2, "tox"]) >= 0.3)
    
    # Initialize selection as NA
    output[Iter, "select dose"] <- NA
    
    # If both are dropped, end trial (no dose selected)
    if (drop_dose1 & drop_dose2) {
      select <- 0
    } else if (drop_dose1 & !drop_dose2) {
      select <- 2
    } else if (!drop_dose1 & drop_dose2) {
      select <- 1
    } else {
      ### Dose selection (if both doses are still viable)
      # Calculate differences
      benefit_diff <- mean(data_S2[data_S2$Group == 2, "orr"]) - mean(data_S2[data_S2$Group == 1, "orr"])
      risk_diff <- mean(data_S2[data_S2$Group == 2, "tox"]) - mean(data_S2[data_S2$Group == 1, "tox"])
      # BR Score = Benefit - Risk
      BR_score <- benefit_diff - risk_diff
      select <- ifelse(BR_score > 0, 2, 1) # if BR_score > 0, then 2 == high dose; otherwise 1 == low dose
    }
    
    output[Iter,"select dose"] = select
    
    
    ###drop one arm and reset the selected group as Group==1
    data_s = subset(data, (Group == 0 | Group == select), select = -c(orr, tox))
    data_s$Group[which(data_s$Group==select)]=1
    if(select==0){next}
    
    ####cut/censor the survival data by events information rate
    ##need to consider the time difference between the enrollment of N1 (Stage1) and N2 (Stage2)
    ##Stage1 subjects have longer time on study at the analysis
    #create a pseudo time reference (calender time from the first enrollment); 
    data_s$dif=NA # the difference between subject enroll time from the first enrollment (t=0)
    data_s$dif[which(data_s$Stage==1)] = 6+mfup-data_s$flup1[which(data_s$Stage==1)]; 
    #---so firstly 10 per arm until dose selection (for fup+D month); 
    #then enroll 20 subject per arm until reach N subjects per arm
    ## after enrolling the last subject for stage 1, there's a wait for fup+D month for dose selection DBL, 
    ## during which the enrollment is paused. 
    remd = N2%%20; quot= N2%/%20#remainder and the quotient :
    data_s$dif[which(data_s$Group==0& data_s$Stage==2)]= sample(c(13+rep(1:quot-1,each=20),+rep(13+quot,remd)),N2)
    data_s$dif[which(data_s$Group==1& data_s$Stage==2)]= sample(c(13+rep(1:quot-1,each=20),+rep(13+quot,remd)),N2)
    
    data_s$ref.OS = data_s$Time+ data_s$dif #calender time from the first enrollment in phase 2/3
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
    dat1 = sort[c("Time","look.OS","Group","Stage","dif")]; 
    colnames(dat1) =c("Time","look","Group","Stage","dif")
    dat1$ind= ifelse(dat1$look %in% 1 , 1, 0 ); sum(dat1$ind)
    dat1$flup=tcut1-dat1$dif; 
    dat1$Time[which(dat1$ind==0)] = dat1$flup[which(dat1$ind==0)];
    dat1 <- dat1[dat1$Time >= 0, ]
    pval1.OS = surv.test(dat1, altern = "greater") ## Superiority
    
    dat1.OS = dat1
    ##OS test for Stage 2 subjects
    pval1_S2_dat = dat1.OS[which(dat1.OS$Stage==2),]
    pval1_S2 = surv.test(pval1_S2_dat, altern = "greater")
    ##OS test for Stage 1 subjects
    pval1_S1_dat = dat1.OS[which(dat1.OS$Stage==1),]
    pval1_S1 = surv.test(pval1_S1_dat, altern = "greater") 
    ######Adjustment and pvalue combination
    pval1_S1_adjD =dunnett(pval1_S1, m=2)
    ###Use sample size N as weights
    pval1_fD=wInvnorm(c(pval1_S1_adjD,pval1_S2),c(nrow(pval1_S1_dat),nrow(pval1_S2_dat)))
    
    #############################################################################################
    ###################
    ##at OS final look (100% OS events)
    dat2 = sort[c("Time","look.OS","Group","Stage","dif")]; 
    colnames(dat2) =c("Time","look","Group","Stage", "dif")
    dat2$ind= ifelse(dat2$look %in% 1:2 , 1, 0 ); sum(dat2$ind)
    dat2$flup=tcut2-dat2$dif; head(dat2);##add censored values for censored case
    dat2$Time[which(dat2$ind==0)] = dat2$flup[which(dat2$ind==0)];
    pval2.OS = surv.test(dat2, altern = "greater") ## Superiority
    
    dat2.OS = dat2;
    ##OS test for Stage 2 subjects
    pval2_S2_dat = dat2.OS[which(dat2.OS$Stage==2),]
    pval2_S2 = surv.test(pval2_S2_dat, altern = "greater")
    ##OS test for Stage 1 subjects
    pval2_S1_dat = dat2.OS[which(dat2.OS$Stage==1),]
    pval2_S1 = surv.test(pval2_S1_dat, altern = "greater") 
    ######Adjustment and pvalue combination
    pval2_S1_adjD =dunnett(pval2_S1, m=2)
    ###Use sample size N as weights
    pval2_fD=wInvnorm(c(pval2_S1_adjD,pval2_S2),c(nrow(pval2_S1_dat),nrow(pval2_S2_dat)))
    
    ##SAVE
    Allps_noadj= c(pval1.OS,pval2.OS)
    Allps_fD= c(pval1_fD,pval2_fD)
    output[Iter,c("IA NoAdj")]= as.numeric(Allps_noadj[1]<pcut.OS[1])
    output[Iter,c("IA Dunnett")]= as.numeric(Allps_fD[1]<pcut.OS[1])
    output[Iter,c("NoAdj")] =as.numeric(sum(Allps_noadj<pcut.OS)>0)
    output[Iter,c("Dunnett")] =as.numeric(sum(Allps_fD<pcut.OS)>0)
    output[Iter,c("decision NoAdj")]= output[Iter,c("NoAdj")]*expanse*(select>0)
    output[Iter,c("decision Dunnett")]= output[Iter,c("Dunnett")]*expanse*(select>0)
  } #end of M
  output_list[[HH]] <- output
} #end of HH

# Dunnett overall power
mean(output_list[["H1"]][["decision Dunnett"]])
# Dunnett overall type I error
mean(output_list[["H0"]][["decision Dunnett"]])
# Non adjusted overall power
mean(output_list[["H1"]][["decision NoAdj"]])
# Non adjusted overall type I error
mean(output_list[["H0"]][["decision NoAdj"]])

# Dunnett phase 2/3 power
sum(output_list[["H1"]][["decision Dunnett"]])/sum(output_list[["H1"]][["expansion"]])
# Dunnett phase 2/3 type I error
sum(output_list[["H0"]][["decision Dunnett"]])/sum(output_list[["H0"]][["expansion"]])
# Non adjusted phase 2/3 power
sum(output_list[["H1"]][["decision NoAdj"]])/sum(output_list[["H1"]][["expansion"]])
# Non adjusted phase 2/3 type I error
sum(output_list[["H0"]][["decision NoAdj"]])/sum(output_list[["H0"]][["expansion"]])



output_suc = output[output$`decision Dunnett`!=0,]
output_suc <- output_suc %>%
  mutate(
    ess = case_when(
      `IA Dunnett` == 0 ~ 890,
      `IA Dunnett` == 1 & (`IA time` + 2) >= 30 ~ 890,
      `IA Dunnett` == 1 & (`IA time` + 2) >= 29 ~ 870,
      `IA Dunnett` == 1 & (`IA time` + 2) >= 28 ~ 830,
      `IA Dunnett` == 1 & (`IA time` + 2) >= 27 ~ 790,
      `IA Dunnett` == 1 & (`IA time` + 2) >= 26 ~ 750,
      `IA Dunnett` == 1 & (`IA time` + 2) >= 25 ~ 710,
      `IA Dunnett` == 1 & (`IA time` + 2) >= 24 ~ 670,
      TRUE ~ NA_real_
    )
  )
ess1 = mean(output_suc$ess)
duration = sum(output_suc$`FA time`*(1-output_suc$`IA Dunnett`)+output_suc$`IA time`*output_suc$`IA Dunnett`)/dim(output_suc)[1]
duration1 = duration+2 # read out 2 months
duration1
mean(output_suc$`FA time`) # then 2 months read out follow it
mean(output_suc$`IA time`) # then 2 months read out follow it

output_fail = output[output$`decision Dunnett`==0,]
output_fail$ess[output_fail$`DS time`==6] = 20*2
output_fail$ess[output_fail$`DS time`==13&output_fail$`FA time` == 0] = 20*2+50*3 
output_fail$ess[output_fail$`FA time` != 0] = 890 
ess0 = mean(output_fail$ess)
idx <- which(output_fail$`FA time` == 0)
output_fail$`FA time`[idx] <- output_fail$`DS time`[idx]
output_fail$`FA time`[-idx] <- output_fail$`FA time`[-idx]+2 # read out 2 months
duration0 = mean(output_fail$`FA time`)
duration0

table(output_list[["H1"]][["expansion"]])/M
table(output_list[["H0"]][["expansion"]])/M

table(output[["select dose"]])/M
table(output[["select dose"]])/M

table(output[["select dose"]],output[["expansion"]])/M


# when HR_low = 1, HR_high = 0.7
# Dunnett overall power
mean(output[["decision Dunnett"]]&output[["select dose"]]==2)
# Dunnett phase 2/3 power
sum(output[["decision Dunnett"]]==1&output[["select dose"]]==2)/sum(output[["expansion"]])


HR_string <- paste(HR, collapse = "_")
save.image(file = paste0("S02_N2_",N2,"_HR_",HR_string,".RData"))